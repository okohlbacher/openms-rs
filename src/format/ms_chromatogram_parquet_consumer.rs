// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writes OpenSWATH chromatograms to a Parquet file (`.xic`).
//!
//! Ports `FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h` /
//! `MSChromatogramParquetConsumer.cpp`. Each consumed chromatogram becomes
//! one row of the PyProphet-compatible `.xic` table
//! ([`XICSchema`](crate::format::arrow_schema_registry::XICSchema)), annotated
//! from the transition library the writer was built with, and
//! [`XicParquetFile`](crate::format::xic_parquet_file::XicParquetFile) reads
//! it back.
//!
//! | Column | Type | Content |
//! |---|---|---|
//! | `RUN_ID` | int64 | the run id given to [`MSChromatogramParquetConsumer::new`] |
//! | `SOURCE_FILE` | string | the source file given there (null when empty) |
//! | `MS_LEVEL` | int64 | 1 for a precursor trace, 2 for a fragment trace |
//! | `PRECURSOR_ID` | int64 | the compound's numeric id |
//! | `TRANSITION_ID` | int64, null for precursors | the transition's numeric id |
//! | `MODIFIED_SEQUENCE` | string | the compound's sequence |
//! | `PRECURSOR_CHARGE` | int64, null unless positive | the compound's charge |
//! | `PRODUCT_CHARGE` | int64, null unless positive | the fragment charge |
//! | `DETECTING_TRANSITION` | int64 | 1 or 0 |
//! | `PRECURSOR_DECOY` | int64 | 1 when the precursor has only decoy transitions |
//! | `PRODUCT_DECOY` | int64 | the transition's decoy flag |
//! | `TRANSITION_ORDINAL` | int64, null when negative | the fragment number |
//! | `TRANSITION_TYPE` | string | the fragment type (`y`, `b`, ...) |
//! | `ANNOTATION` | string | `y3^1`, or `Precursor_i0` for a precursor trace |
//! | `RT_DATA`, `INTENSITY_DATA` | binary | the encoded arrays |
//! | `RT_COMPRESSION`, `INTENSITY_COMPRESSION` | int64 | 5 and 6 (see below) |
//!
//! Numeric ids: a compound id or transition name that is an integer is its
//! own id; any other is numbered from 1 upwards past the ids already taken.
//!
//! Arrays: retention times are MS-Numpress linear (fixed point for a mass
//! accuracy of 0.05) and intensities MS-Numpress slof, each then
//! zlib-compressed, identifiers 5 and 6, as in the source. Where the lossy
//! encoding would not decode back to the data, the array is stored as
//! zlib-compressed doubles under identifier 1 instead (native; see
//! `docs/OPENSWATH_PARQUET_SUPPORT.md`): in particular a chromatogram of one
//! or two points, which the source's encoder stores with a fixed point of 0.
//!
//! ```no_run
//! use openms::analysis::openswath::transition_experiment::LightTargetedExperiment;
//! use openms::format::ms_chromatogram_parquet_consumer::MSChromatogramParquetConsumer;
//! use openms::MSChromatogram;
//!
//! # fn main() -> openms::Result<()> {
//! let library = LightTargetedExperiment::default();
//! let mut consumer = MSChromatogramParquetConsumer::new("out.xic", 1, "run.mzML", &library);
//! let mut chromatogram = MSChromatogram::default();
//! chromatogram.native_id = "tr1".into();
//! # let _ = &mut chromatogram;
//! // consumer.consume(&mut chromatogram)?; // needs "tr1" in the library
//! consumer.finalize()?;
//! # Ok(())
//! # }
//! ```

use crate::analysis::openswath::openswath_helper::compute_transition_group_id;
use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::format::arrow_schema_registry::XICSchema;
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::numpress_coder::NumpressCoderLimits;
use crate::format::openswath_parquet_writer::{
    self as writer, Lossy, StagedParquet, TransitionMaps,
};
use crate::interfaces::MSDataConsumer;
use crate::metadata::ExperimentalSettings;
use crate::{Error, MSChromatogram, MSSpectrum, Result};
use arrow_array::builder::{BinaryBuilder, Int64Builder, StringBuilder};
use arrow_array::{ArrayRef, RecordBatch};
use std::ops::ControlFlow;
use std::path::Path;
use std::sync::Arc;

/// The source's `flush_row_limit_`: rows buffered before a row group is
/// written.
pub const FLUSH_ROW_LIMIT: usize = 1000;
/// The source's `flush_byte_limit_`: encoded bytes buffered before a row
/// group is written (below Arrow's 2 GiB binary-array limit).
pub const FLUSH_BYTE_LIMIT: usize = 1_900_000_000;
/// Row group size of the source's `WriteTable(table, 1024)`.
const ROW_GROUP_ROWS: usize = 1024;
/// Mass accuracy the source's retention-time encoder uses.
const RT_MASS_ACCURACY: f64 = 0.05;
const PRECURSOR_TAG: &str = "_Precursor_i";

/// The metadata part of one row, resolved from the library.
#[derive(Clone, Debug, Default)]
struct Row {
    ms_level: i64,
    precursor_id: Option<i64>,
    transition_id: Option<i64>,
    modified_sequence: String,
    precursor_charge: Option<i64>,
    product_charge: Option<i64>,
    detecting_transition: Option<i64>,
    precursor_decoy: Option<i64>,
    product_decoy: Option<i64>,
    transition_ordinal: Option<i64>,
    transition_type: String,
    annotation: String,
    rt: writer::Encoded,
    intensity: writer::Encoded,
}

/// Writes chromatograms to a Parquet file with a PyProphet-compatible schema.
///
/// Ports `OpenMS::MSChromatogramParquetConsumer`, an `IMSDataConsumer`:
/// spectra are ignored (their peaks are cleared), each chromatogram is
/// encoded and buffered, and every [`FLUSH_ROW_LIMIT`] rows (or
/// [`FLUSH_BYTE_LIMIT`] encoded bytes) a row group is written.
/// [`finalize`](Self::finalize) writes the rest and closes the file; a
/// consumer that was never given a chromatogram still writes a file with no
/// rows.
///
/// Dropping an unfinalised consumer finalises it and discards any error, as
/// the source's destructor does after logging it; call `finalize` to see
/// errors.
///
/// Native differences: the file is staged next to the destination and moved
/// onto it by `finalize`, so the destination is untouched until then and a
/// failure never leaves a partial file (the source truncates it at the first
/// flush); and a consumer refuses chromatograms after `finalize` (the
/// source's would buffer them and, 1000 rows later, reopen and truncate the
/// finished file).
pub struct MSChromatogramParquetConsumer {
    run_id: i64,
    source_file: String,
    maps: TransitionMaps,
    pending: Vec<Row>,
    pending_bytes: usize,
    output: Option<StagedParquet>,
    finished: bool,
    /// Ceilings for encoding one array (native).
    pub limits: NumpressCoderLimits,
}

impl std::fmt::Debug for MSChromatogramParquetConsumer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MSChromatogramParquetConsumer")
            .field("run_id", &self.run_id)
            .field("source_file", &self.source_file)
            .field("pending", &self.pending.len())
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl MSChromatogramParquetConsumer {
    /// A consumer writing to `filename`.
    ///
    /// `run_id` is stored in every row's `RUN_ID`; a value above `i64::MAX`
    /// wraps to a negative one, as the source's cast does. The library
    /// annotates the chromatograms: a fragment trace is looked up by its
    /// native id among the transition names, a precursor trace (native id
    /// containing `_Precursor_i`) by the transition group id of its native id
    /// among the compounds. Nothing is written before the first flush.
    pub fn new(
        filename: impl AsRef<Path>,
        run_id: u64,
        source_file: &str,
        transition_exp: &LightTargetedExperiment,
    ) -> Self {
        Self {
            run_id: run_id as i64,
            source_file: source_file.to_owned(),
            maps: TransitionMaps::build(transition_exp, false),
            pending: Vec::new(),
            pending_bytes: 0,
            output: Some(StagedParquet::new(filename.as_ref(), XICSchema::schema())),
            finished: false,
            limits: NumpressCoderLimits::default(),
        }
    }

    /// Resolve and encode one chromatogram, then clear its peaks (the
    /// source's `consumeChromatogram`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the native id is not in the library
    /// ("does not have a matching transition metadata entry", or, for a
    /// precursor trace, "looks like a precursor but no matching compound
    /// metadata entry was found"), after `finalize`, and when a flush fails;
    /// [`Error::Io`] when the output cannot be created. On an error the
    /// chromatogram is left unchanged and nothing is buffered.
    pub fn consume(&mut self, chromatogram: &mut MSChromatogram) -> Result<()> {
        if self.finished {
            return Err(Error::InvalidValue(
                "MSChromatogramParquetConsumer cannot accept chromatograms after finalize()".into(),
            ));
        }
        let native_id = chromatogram.native_id.as_str();
        let mut row = if native_id.contains(PRECURSOR_TAG) {
            let group_id = compute_transition_group_id(native_id);
            let compound = self.maps.compounds.get(&group_id).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Chromatogram native ID '{native_id}' looks like a precursor but no matching \
                     compound metadata entry was found. Please ensure the native ID is correct \
                     and that the transition library contains metadata for this precursor."
                ))
            })?;
            // `buildPrecursorAnnotation_`: the text after the last
            // "_Precursor_i" tag's leading underscore.
            let annotation = native_id
                .rfind(PRECURSOR_TAG)
                .map(|position| native_id[position + 1..].to_owned())
                .unwrap_or_default();
            Row {
                ms_level: 1,
                precursor_id: Some(compound.precursor_id),
                modified_sequence: compound.modified_sequence.clone(),
                precursor_charge: (compound.precursor_charge > 0)
                    .then_some(compound.precursor_charge),
                precursor_decoy: Some(compound.precursor_decoy),
                annotation,
                ..Row::default()
            }
        } else {
            let transition = self.maps.transitions.get(native_id).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Chromatogram native ID '{native_id}' does not have a matching transition \
                     metadata entry. Please ensure the native ID is correct and that the \
                     transition library contains metadata for this transition."
                ))
            })?;
            Row {
                ms_level: 2,
                precursor_id: Some(transition.precursor_id),
                transition_id: Some(transition.transition_id),
                modified_sequence: transition.modified_sequence.clone(),
                precursor_charge: (transition.precursor_charge > 0)
                    .then_some(transition.precursor_charge),
                product_charge: (transition.product_charge > 0)
                    .then_some(transition.product_charge),
                detecting_transition: Some(transition.detecting_transition),
                precursor_decoy: Some(transition.precursor_decoy),
                product_decoy: Some(transition.product_decoy),
                transition_ordinal: (transition.transition_ordinal >= 0)
                    .then_some(transition.transition_ordinal),
                transition_type: transition.transition_type.clone(),
                annotation: transition.annotation.clone(),
                ..Row::default()
            }
        };
        let rt: Vec<f64> = chromatogram.peaks.iter().map(|peak| peak.rt).collect();
        let intensity: Vec<f64> = chromatogram
            .peaks
            .iter()
            .map(|peak| f64::from(peak.intensity))
            .collect();
        row.rt = writer::encode(
            &rt,
            Lossy::Linear {
                mass_accuracy: Some(RT_MASS_ACCURACY),
            },
            self.limits,
        )?;
        row.intensity = writer::encode(&intensity, Lossy::Slof, self.limits)?;
        if row.rt.bytes.len() > i32::MAX as usize || row.intensity.bytes.len() > i32::MAX as usize {
            return Err(Error::InvalidValue(
                "Binary value exceeds Arrow int32 size limit (2 GB).".into(),
            ));
        }
        self.pending_bytes = self
            .pending_bytes
            .saturating_add(row.rt.bytes.len() + row.intensity.bytes.len());
        self.pending.push(row);
        chromatogram.clear(false);
        if self.pending.len() >= FLUSH_ROW_LIMIT || self.pending_bytes >= FLUSH_BYTE_LIMIT {
            self.flush()?;
        }
        Ok(())
    }

    /// Write the buffered rows as row groups (`flushBuffers_`).
    fn flush(&mut self) -> Result<()> {
        let rows = std::mem::take(&mut self.pending);
        self.pending_bytes = 0;
        let batch = self.batch(&rows)?;
        let output = self
            .output
            .as_mut()
            .ok_or_else(|| Error::InvalidValue("parquet writer unavailable".into()))?;
        output.write(&batch, ROW_GROUP_ROWS)
    }

    /// The buffered rows as a table of the XIC schema.
    fn batch(&self, rows: &[Row]) -> Result<RecordBatch> {
        let mut run_id = Int64Builder::new();
        let mut source_file = StringBuilder::new();
        let mut ms_level = Int64Builder::new();
        let mut precursor_id = Int64Builder::new();
        let mut transition_id = Int64Builder::new();
        let mut modified_sequence = StringBuilder::new();
        let mut precursor_charge = Int64Builder::new();
        let mut product_charge = Int64Builder::new();
        let mut detecting_transition = Int64Builder::new();
        let mut precursor_decoy = Int64Builder::new();
        let mut product_decoy = Int64Builder::new();
        let mut transition_ordinal = Int64Builder::new();
        let mut transition_type = StringBuilder::new();
        let mut annotation = StringBuilder::new();
        let mut rt_data = BinaryBuilder::new();
        let mut intensity_data = BinaryBuilder::new();
        let mut rt_compression = Int64Builder::new();
        let mut intensity_compression = Int64Builder::new();
        let text = |builder: &mut StringBuilder, value: &str| {
            builder.append_option((!value.is_empty()).then_some(value));
        };
        for row in rows {
            run_id.append_value(self.run_id);
            text(&mut source_file, &self.source_file);
            ms_level.append_value(row.ms_level);
            precursor_id.append_option(row.precursor_id);
            transition_id.append_option(row.transition_id);
            text(&mut modified_sequence, &row.modified_sequence);
            precursor_charge.append_option(row.precursor_charge);
            product_charge.append_option(row.product_charge);
            detecting_transition.append_option(row.detecting_transition);
            precursor_decoy.append_option(row.precursor_decoy);
            product_decoy.append_option(row.product_decoy);
            transition_ordinal.append_option(row.transition_ordinal);
            text(&mut transition_type, &row.transition_type);
            text(&mut annotation, &row.annotation);
            rt_data.append_value(&row.rt.bytes);
            intensity_data.append_value(&row.intensity.bytes);
            rt_compression.append_value(row.rt.code);
            intensity_compression.append_value(row.intensity.code);
        }
        let columns: Vec<ArrayRef> = vec![
            Arc::new(run_id.finish()),
            Arc::new(source_file.finish()),
            Arc::new(ms_level.finish()),
            Arc::new(precursor_id.finish()),
            Arc::new(transition_id.finish()),
            Arc::new(modified_sequence.finish()),
            Arc::new(precursor_charge.finish()),
            Arc::new(product_charge.finish()),
            Arc::new(detecting_transition.finish()),
            Arc::new(precursor_decoy.finish()),
            Arc::new(product_decoy.finish()),
            Arc::new(transition_ordinal.finish()),
            Arc::new(transition_type.finish()),
            Arc::new(annotation.finish()),
            Arc::new(rt_data.finish()),
            Arc::new(intensity_data.finish()),
            Arc::new(rt_compression.finish()),
            Arc::new(intensity_compression.finish()),
        ];
        let table = RecordBatch::try_new(XICSchema::schema(), columns)
            .map_err(|error| Error::InvalidValue(format!("XIC table: {error}")))?;
        let check = validation::validate(&table, &XICSchema::schema(), Mode::Strict);
        if !check.valid {
            return Err(Error::InvalidValue(format!(
                "XIC table schema validation failed: {check}"
            )));
        }
        Ok(table)
    }

    /// Write the buffered rows and close the file (`finalize`).
    ///
    /// A second call does nothing. A consumer that received no chromatogram
    /// writes a file with the XIC schema and no rows.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the output cannot be created or moved onto the
    /// destination (the source's `Exception::FileNotWritable`);
    /// [`Error::InvalidValue`] when Parquet encoding fails. The consumer is
    /// finished either way.
    pub fn finalize(&mut self) -> Result<()> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        if !self.pending.is_empty() {
            self.flush()?;
        }
        match self.output.take() {
            Some(output) => output.finish(),
            None => Ok(()),
        }
    }
}

impl Drop for MSChromatogramParquetConsumer {
    fn drop(&mut self) {
        let _ = self.finalize();
    }
}

impl MSDataConsumer for MSChromatogramParquetConsumer {
    /// Reserves buffer space for up to [`FLUSH_ROW_LIMIT`] chromatograms;
    /// spectra are ignored.
    fn set_expected_size(&mut self, _spectra: usize, chromatograms: usize) -> Result<()> {
        self.pending
            .try_reserve(chromatograms.min(FLUSH_ROW_LIMIT))
            .map_err(|_| Error::InvalidValue("Failed to reserve capacity".into()))
    }

    /// Unused, as in the source.
    fn set_experimental_settings(&mut self, _settings: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }

    /// Spectra are not exported: the peaks are cleared, the metadata kept.
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        spectrum.clear(false);
        Ok(ControlFlow::Continue(()))
    }

    /// [`consume`](MSChromatogramParquetConsumer::consume).
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        self.consume(chromatogram)?;
        Ok(ControlFlow::Continue(()))
    }
}
