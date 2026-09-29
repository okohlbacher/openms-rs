// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writes OpenSWATH mobilograms (ion-mobility traces) to a Parquet file
//! (`.xim`).
//!
//! Ports `FORMAT/DATAACCESS/MobilogramParquetConsumer.h` /
//! `MobilogramParquetConsumer.cpp`. Each consumed mobilogram becomes one row
//! of the `.xim` table
//! ([`XIMSchema`](crate::format::arrow_schema_registry::XIMSchema)), which
//! [`XimParquetFile`](crate::format::xim_parquet_file::XimParquetFile) reads
//! back. The row carries the run, the caller's annotation
//! ([`MobilogramAnnotation`]) and, when the annotation names a transition of
//! the library, that transition's precursor and fragment metadata, with the
//! same numeric ids and decoy rule as the `.xic` writer
//! ([`ms_chromatogram_parquet_consumer`](crate::format::ms_chromatogram_parquet_consumer)).
//!
//! Arrays: mobilities are MS-Numpress linear (maximal fixed point) and
//! intensities MS-Numpress slof, each then zlib-compressed, identifiers 5 and
//! 6, as in the source; an array the lossy encoding would not decode back is
//! stored as zlib-compressed doubles under identifier 1 (native; see
//! `docs/OPENSWATH_PARQUET_SUPPORT.md`).
//!
//! ```no_run
//! use openms::analysis::openswath::transition_experiment::LightTargetedExperiment;
//! use openms::format::mobilogram_parquet_consumer::{MobilogramAnnotation, MobilogramParquetConsumer};
//! use openms::{MobilityPeak1D, Mobilogram};
//!
//! # fn main() -> openms::Result<()> {
//! let mut consumer =
//!     MobilogramParquetConsumer::new("out.xim", 1, "run.mzML", &LightTargetedExperiment::default());
//! let mobilogram = Mobilogram::from_peaks(vec![
//!     MobilityPeak1D::new(0.8, 1000.0),
//!     MobilityPeak1D::new(0.9, 2000.0),
//! ]);
//! let annotation = MobilogramAnnotation {
//!     mobilogram_type: "ms1".into(),
//!     ms_level: Some(1),
//!     ..MobilogramAnnotation::default()
//! };
//! consumer.consume_mobilogram(&mobilogram, &annotation)?;
//! consumer.finalize()?;
//! # Ok(())
//! # }
//! ```

use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::format::arrow_schema_registry::XIMSchema;
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::numpress_coder::NumpressCoderLimits;
use crate::format::openswath_parquet_writer::{
    self as writer, Lossy, StagedParquet, TransitionMaps,
};
use crate::kernel::Mobilogram;
use crate::{Error, Result};
use arrow_array::builder::{BinaryBuilder, Float64Builder, Int64Builder, StringBuilder};
use arrow_array::{ArrayRef, RecordBatch};
use std::path::Path;
use std::sync::Arc;

/// The source's `MAX_BUFFERED_BINARY_BYTES_`: encoded bytes buffered before
/// the rows are written.
pub const MAX_BUFFERED_BINARY_BYTES: usize = 256 * 1024 * 1024;
/// The source's `PARQUET_ROW_GROUP_SIZE_`.
const ROW_GROUP_ROWS: usize = 1024;

/// What the caller says about one mobilogram: the optional arguments of the
/// source's `consumeMobilogram`.
///
/// The source's sentinels are `""`, `-1` and NaN; here `""`, `None` and
/// `None` (a negative integer or a non-finite RT is unset too, as in the
/// source).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MobilogramAnnotation {
    /// `MOBILOGRAM_TYPE`, for example `ms1`, `ms2` or `ms1_contrast`.
    pub mobilogram_type: String,
    /// `MS_LEVEL`.
    pub ms_level: Option<i64>,
    /// `TRANSITION_ID`; unset, the library transition's id is used.
    pub transition_id: Option<i64>,
    /// The native id of the library transition this mobilogram belongs to;
    /// its precursor and fragment metadata fill the row.
    pub transition_native_id: String,
    /// `FEATURE_RT`, the RT apex the mobilogram belongs to.
    pub feature_rt: Option<f64>,
    /// `FEATURE_ID`.
    pub feature_id: Option<i64>,
}

#[derive(Clone, Debug, Default)]
struct Row {
    ms_level: Option<i64>,
    mobilogram_type: String,
    precursor_id: Option<i64>,
    transition_id: Option<i64>,
    feature_id: Option<i64>,
    feature_rt: Option<f64>,
    modified_sequence: String,
    precursor_charge: Option<i64>,
    product_charge: Option<i64>,
    detecting_transition: Option<i64>,
    precursor_decoy: Option<i64>,
    product_decoy: Option<i64>,
    transition_ordinal: Option<i64>,
    transition_type: String,
    annotation: String,
    mobility: writer::Encoded,
    intensity: writer::Encoded,
}

/// Writes mobilograms to a Parquet file.
///
/// Ports `OpenMS::MobilogramParquetConsumer`. Rows are buffered until their
/// encoded arrays reach [`MAX_BUFFERED_BINARY_BYTES`], then written as row
/// groups of 1024; [`finalize`](Self::finalize) writes the rest and closes
/// the file, and a consumer that was never given a mobilogram still writes a
/// file with no rows. The source serialises concurrent `consumeMobilogram`
/// calls with a mutex; here the writer takes `&mut self`, so a caller that
/// shares one across threads wraps it in a `Mutex`.
///
/// Dropping an unfinalised consumer finalises it and discards any error, as
/// the source's destructor does after logging it; call `finalize` to see
/// errors. The file is staged next to the destination and moved onto it by
/// `finalize` (native: the source truncates the destination at its first
/// write).
pub struct MobilogramParquetConsumer {
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

impl std::fmt::Debug for MobilogramParquetConsumer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MobilogramParquetConsumer")
            .field("run_id", &self.run_id)
            .field("source_file", &self.source_file)
            .field("pending", &self.pending.len())
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl MobilogramParquetConsumer {
    /// A consumer writing to `filename`.
    ///
    /// `run_id` is stored with its top bit cleared (`run_id & !(1 << 63)`),
    /// as in the source, so it is never negative. Pass an empty library for
    /// no transition annotation. Nothing is written before the first flush.
    pub fn new(
        filename: impl AsRef<Path>,
        run_id: u64,
        source_file: &str,
        transition_exp: &LightTargetedExperiment,
    ) -> Self {
        Self {
            run_id: (run_id & !(1u64 << 63)) as i64,
            source_file: source_file.to_owned(),
            maps: TransitionMaps::build(transition_exp, true),
            pending: Vec::new(),
            pending_bytes: 0,
            output: Some(StagedParquet::new(filename.as_ref(), XIMSchema::schema())),
            finished: false,
            limits: NumpressCoderLimits::default(),
        }
    }

    /// Reserve buffer space (`setExpectedSize`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the allocation fails.
    pub fn set_expected_size(&mut self, expected_mobilograms: usize) -> Result<()> {
        self.pending
            .try_reserve(expected_mobilograms)
            .map_err(|_| Error::InvalidValue("Failed to reserve capacity".into()))
    }

    /// Encode one mobilogram and buffer it (`consumeMobilogram`). The
    /// mobilogram is not modified.
    ///
    /// When `annotation.transition_native_id` names a library transition, the
    /// row gets its precursor id, sequence, charges, decoy flags, ordinal,
    /// type and annotation, and its transition id unless the annotation gives
    /// one; otherwise those columns are null.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] after `finalize` ("cannot accept writes after
    /// finalize()", the source's `Exception::IllegalArgument`) and when a
    /// flush fails; [`Error::Io`] when the output cannot be created. On an
    /// error nothing is buffered.
    pub fn consume_mobilogram(
        &mut self,
        mobilogram: &Mobilogram,
        annotation: &MobilogramAnnotation,
    ) -> Result<()> {
        if self.finished {
            return Err(Error::InvalidValue(
                "MobilogramParquetConsumer cannot accept writes after finalize()".into(),
            ));
        }
        let mobility: Vec<f64> = mobilogram.peaks.iter().map(|peak| peak.mobility).collect();
        let intensity: Vec<f64> = mobilogram
            .peaks
            .iter()
            .map(|peak| f64::from(peak.intensity))
            .collect();
        let mobility = writer::encode(
            &mobility,
            Lossy::Linear {
                mass_accuracy: None,
            },
            self.limits,
        )?;
        let intensity = writer::encode(&intensity, Lossy::Slof, self.limits)?;
        if mobility.bytes.len() > i32::MAX as usize || intensity.bytes.len() > i32::MAX as usize {
            return Err(Error::InvalidValue(
                "Single binary payload exceeds Arrow binary cell limit".into(),
            ));
        }
        let size = mobility.bytes.len() + intensity.bytes.len();
        if !self.pending.is_empty()
            && self.pending_bytes.saturating_add(size) > MAX_BUFFERED_BINARY_BYTES
        {
            self.flush()?;
        }
        let nonnegative = |value: Option<i64>| value.filter(|value| *value >= 0);
        let mut row = Row {
            ms_level: nonnegative(annotation.ms_level),
            mobilogram_type: annotation.mobilogram_type.clone(),
            feature_id: nonnegative(annotation.feature_id),
            feature_rt: annotation.feature_rt.filter(|rt| rt.is_finite()),
            mobility,
            intensity,
            ..Row::default()
        };
        match self.maps.transitions.get(&annotation.transition_native_id) {
            Some(transition) => {
                row.precursor_id = Some(transition.precursor_id);
                row.transition_id =
                    Some(nonnegative(annotation.transition_id).unwrap_or(transition.transition_id));
                row.modified_sequence = transition.modified_sequence.clone();
                row.precursor_charge =
                    (transition.precursor_charge > 0).then_some(transition.precursor_charge);
                row.product_charge =
                    (transition.product_charge > 0).then_some(transition.product_charge);
                row.detecting_transition = Some(transition.detecting_transition);
                row.precursor_decoy = Some(transition.precursor_decoy);
                row.product_decoy = Some(transition.product_decoy);
                row.transition_ordinal =
                    (transition.transition_ordinal >= 0).then_some(transition.transition_ordinal);
                row.transition_type = transition.transition_type.clone();
                row.annotation = transition.annotation.clone();
            }
            None => row.transition_id = nonnegative(annotation.transition_id),
        }
        self.pending_bytes = self.pending_bytes.saturating_add(size);
        self.pending.push(row);
        if self.pending_bytes >= MAX_BUFFERED_BINARY_BYTES {
            self.flush()?;
        }
        Ok(())
    }

    /// `flushChunk_`.
    fn flush(&mut self) -> Result<()> {
        if self.pending.is_empty() {
            return Ok(());
        }
        let rows = std::mem::take(&mut self.pending);
        self.pending_bytes = 0;
        let batch = self.batch(&rows)?;
        let output = self
            .output
            .as_mut()
            .ok_or_else(|| Error::InvalidValue("parquet writer unavailable".into()))?;
        output.write(&batch, ROW_GROUP_ROWS)
    }

    /// `buildTableFromBuilders_`.
    fn batch(&self, rows: &[Row]) -> Result<RecordBatch> {
        let mut run_id = Int64Builder::new();
        let mut source_file = StringBuilder::new();
        let mut ms_level = Int64Builder::new();
        let mut mobilogram_type = StringBuilder::new();
        let mut precursor_id = Int64Builder::new();
        let mut transition_id = Int64Builder::new();
        let mut feature_id = Int64Builder::new();
        let mut feature_rt = Float64Builder::new();
        let mut modified_sequence = StringBuilder::new();
        let mut precursor_charge = Int64Builder::new();
        let mut product_charge = Int64Builder::new();
        let mut detecting_transition = Int64Builder::new();
        let mut precursor_decoy = Int64Builder::new();
        let mut product_decoy = Int64Builder::new();
        let mut transition_ordinal = Int64Builder::new();
        let mut transition_type = StringBuilder::new();
        let mut annotation = StringBuilder::new();
        let mut mobility_data = BinaryBuilder::new();
        let mut intensity_data = BinaryBuilder::new();
        let mut mobility_compression = Int64Builder::new();
        let mut intensity_compression = Int64Builder::new();
        let text = |builder: &mut StringBuilder, value: &str| {
            builder.append_option((!value.is_empty()).then_some(value));
        };
        for row in rows {
            run_id.append_value(self.run_id);
            text(&mut source_file, &self.source_file);
            ms_level.append_option(row.ms_level);
            text(&mut mobilogram_type, &row.mobilogram_type);
            precursor_id.append_option(row.precursor_id);
            transition_id.append_option(row.transition_id);
            feature_id.append_option(row.feature_id);
            feature_rt.append_option(row.feature_rt);
            text(&mut modified_sequence, &row.modified_sequence);
            precursor_charge.append_option(row.precursor_charge);
            product_charge.append_option(row.product_charge);
            detecting_transition.append_option(row.detecting_transition);
            precursor_decoy.append_option(row.precursor_decoy);
            product_decoy.append_option(row.product_decoy);
            transition_ordinal.append_option(row.transition_ordinal);
            text(&mut transition_type, &row.transition_type);
            text(&mut annotation, &row.annotation);
            mobility_data.append_value(&row.mobility.bytes);
            intensity_data.append_value(&row.intensity.bytes);
            mobility_compression.append_value(row.mobility.code);
            intensity_compression.append_value(row.intensity.code);
        }
        let columns: Vec<ArrayRef> = vec![
            Arc::new(run_id.finish()),
            Arc::new(source_file.finish()),
            Arc::new(ms_level.finish()),
            Arc::new(mobilogram_type.finish()),
            Arc::new(precursor_id.finish()),
            Arc::new(transition_id.finish()),
            Arc::new(feature_id.finish()),
            Arc::new(feature_rt.finish()),
            Arc::new(modified_sequence.finish()),
            Arc::new(precursor_charge.finish()),
            Arc::new(product_charge.finish()),
            Arc::new(detecting_transition.finish()),
            Arc::new(precursor_decoy.finish()),
            Arc::new(product_decoy.finish()),
            Arc::new(transition_ordinal.finish()),
            Arc::new(transition_type.finish()),
            Arc::new(annotation.finish()),
            Arc::new(mobility_data.finish()),
            Arc::new(intensity_data.finish()),
            Arc::new(mobility_compression.finish()),
            Arc::new(intensity_compression.finish()),
        ];
        let table = RecordBatch::try_new(XIMSchema::schema(), columns)
            .map_err(|error| Error::InvalidValue(format!("XIM table: {error}")))?;
        let check = validation::validate(&table, &XIMSchema::schema(), Mode::Strict);
        if !check.valid {
            return Err(Error::InvalidValue(format!(
                "XIM table schema validation failed: {check}"
            )));
        }
        Ok(table)
    }

    /// Write the buffered rows and close the file (`finalize`).
    ///
    /// A second call does nothing. A consumer that received no mobilogram
    /// writes a file with the XIM schema and no rows.
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
        self.flush()?;
        match self.output.take() {
            Some(output) => output.finish(),
            None => Ok(()),
        }
    }
}

impl Drop for MobilogramParquetConsumer {
    fn drop(&mut self) {
        let _ = self.finalize();
    }
}
