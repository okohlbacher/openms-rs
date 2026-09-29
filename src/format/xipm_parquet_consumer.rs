// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Writes targeted m/z–RT–ion-mobility peak maps to a Parquet file (`.xipm`).
//!
//! Ports `FORMAT/DATAACCESS/XIPMParquetConsumer.h` /
//! `XIPMParquetConsumer.cpp`. Each consumed
//! [`ExtractedPeakMap`](crate::analysis::openswath::peak_map_extractor::ExtractedPeakMap)
//! becomes one row of the `.xipm` table
//! ([`XIPMSchema`](crate::format::arrow_schema_registry::XIPMSchema)): run
//! metadata, the precursor or transition annotation from the library (with
//! the same numeric ids and decoy rule as the `.xic` writer), the extraction
//! target and window, and the four arrays. [`XipmParquetFile`](crate::format::xipm_parquet_file::XipmParquetFile)
//! reads it back.
//!
//! A peak map whose native id contains `_Precursor_i` is a `precursor` row,
//! looked up by its transition group id among the compounds; any other is a
//! `transition` row, looked up among the transition names. `TARGET_RT` is the
//! peak map's target RT, or the midpoint of a non-empty RT window;
//! `RT_START`/`RT_END` are set only for a non-empty window, and
//! `TARGET_ION_MOBILITY` only when not negative.
//!
//! Arrays: m/z, RT and ion mobility are MS-Numpress linear (maximal fixed
//! point), intensities MS-Numpress slof, each then zlib-compressed,
//! identifiers 5, 5, 5 and 6, as in the source; an array the lossy encoding
//! would not decode back is stored as zlib-compressed doubles under
//! identifier 1 (native; see `docs/OPENSWATH_PARQUET_SUPPORT.md`).

use crate::analysis::openswath::openswath_helper::compute_transition_group_id;
use crate::analysis::openswath::peak_map_extractor::ExtractedPeakMap;
use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::format::arrow_schema_registry::XIPMSchema;
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::numpress_coder::NumpressCoderLimits;
use crate::format::openswath_parquet_writer::{
    self as writer, Lossy, StagedParquet, TransitionMaps,
};
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
const PRECURSOR_TAG: &str = "_Precursor_i";
const LINEAR: Lossy = Lossy::Linear {
    mass_accuracy: None,
};

#[derive(Clone, Debug, Default)]
struct Row {
    run_id: i64,
    source_file: String,
    ms_level: i64,
    peakmap_type: &'static str,
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
    target_mz: f64,
    target_rt: Option<f64>,
    target_ion_mobility: Option<f64>,
    rt_start: Option<f64>,
    rt_end: Option<f64>,
    arrays: [writer::Encoded; 4],
}

/// Writes targeted peak maps to a Parquet file.
///
/// Ports `OpenMS::XIPMParquetConsumer`. Rows are buffered until their encoded
/// arrays reach [`MAX_BUFFERED_BINARY_BYTES`], then written as row groups of
/// 1024; [`finalize`](Self::finalize) writes the rest and closes the file,
/// and a consumer that was never given a peak map still writes a file with
/// no rows. Unlike the chromatogram and mobilogram writers, the run id and
/// source file are given per peak map.
///
/// Dropping an unfinalised consumer finalises it and discards any error, as
/// the source's destructor does after logging it; call `finalize` to see
/// errors. The file is staged next to the destination and moved onto it by
/// `finalize` (native: the source truncates the destination at its first
/// write).
pub struct XIPMParquetConsumer {
    maps: TransitionMaps,
    pending: Vec<Row>,
    pending_bytes: usize,
    output: Option<StagedParquet>,
    finished: bool,
    /// Ceilings for encoding one array (native).
    pub limits: NumpressCoderLimits,
}

impl std::fmt::Debug for XIPMParquetConsumer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("XIPMParquetConsumer")
            .field("pending", &self.pending.len())
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl XIPMParquetConsumer {
    /// A consumer writing to `filename`, annotating from `transition_exp`.
    /// Nothing is written before the first flush.
    pub fn new(filename: impl AsRef<Path>, transition_exp: &LightTargetedExperiment) -> Self {
        Self {
            maps: TransitionMaps::build(transition_exp, false),
            pending: Vec::new(),
            pending_bytes: 0,
            output: Some(StagedParquet::new(filename.as_ref(), XIPMSchema::schema())),
            finished: false,
            limits: NumpressCoderLimits::default(),
        }
    }

    /// Reserve buffer space (`setExpectedSize`).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the allocation fails.
    pub fn set_expected_size(&mut self, expected_peak_maps: usize) -> Result<()> {
        self.pending
            .try_reserve(expected_peak_maps)
            .map_err(|_| Error::InvalidValue("Failed to reserve capacity".into()))
    }

    /// Encode one peak map and buffer it (`consumePeakMap`).
    ///
    /// `run_id` is stored as the source's cast stores it (a value above
    /// `i64::MAX` wraps), `ms_level` as given.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the four arrays differ in length ("All
    /// peak-map arrays need to have the same size.", the source's
    /// `Exception::IllegalArgument`), after `finalize`, when the native id is
    /// not in the library ("does not have a matching transition metadata
    /// entry" / "looks like a precursor but no matching compound metadata
    /// entry was found") and when a flush fails; [`Error::Io`] when the
    /// output cannot be created. On an error nothing is buffered.
    pub fn consume_peak_map(
        &mut self,
        peak_map: &ExtractedPeakMap,
        run_id: u64,
        source_file: &str,
        ms_level: i64,
    ) -> Result<()> {
        let n = peak_map.mz.len();
        if peak_map.rt.len() != n
            || peak_map.ion_mobility.len() != n
            || peak_map.intensity.len() != n
        {
            return Err(Error::InvalidValue(
                "All peak-map arrays need to have the same size.".into(),
            ));
        }
        if self.finished {
            return Err(Error::InvalidValue(
                "XIPMParquetConsumer cannot accept writes after finalize().".into(),
            ));
        }
        let native_id = peak_map.native_id.as_str();
        let mut row = if native_id.contains(PRECURSOR_TAG) {
            let group_id = compute_transition_group_id(native_id);
            let compound = self.maps.compounds.get(&group_id).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Peak map native ID '{native_id}' looks like a precursor but no matching \
                     compound metadata entry was found."
                ))
            })?;
            Row {
                peakmap_type: "precursor",
                precursor_id: Some(compound.precursor_id),
                modified_sequence: compound.modified_sequence.clone(),
                precursor_charge: (compound.precursor_charge > 0)
                    .then_some(compound.precursor_charge),
                precursor_decoy: Some(compound.precursor_decoy),
                annotation: native_id
                    .rfind(PRECURSOR_TAG)
                    .map(|position| native_id[position + 1..].to_owned())
                    .unwrap_or_default(),
                ..Row::default()
            }
        } else {
            let transition = self.maps.transitions.get(native_id).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Peak map native ID '{native_id}' does not have a matching transition \
                     metadata entry."
                ))
            })?;
            Row {
                peakmap_type: "transition",
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
        let arrays = [
            writer::encode(&peak_map.mz, LINEAR, self.limits)?,
            writer::encode(&peak_map.rt, LINEAR, self.limits)?,
            writer::encode(&peak_map.ion_mobility, LINEAR, self.limits)?,
            writer::encode(&peak_map.intensity, Lossy::Slof, self.limits)?,
        ];
        if arrays
            .iter()
            .any(|array| array.bytes.len() > i32::MAX as usize)
        {
            return Err(Error::InvalidValue(
                "Single binary payload exceeds Arrow binary cell limit".into(),
            ));
        }
        let size: usize = arrays.iter().map(|array| array.bytes.len()).sum();
        let has_window = peak_map.rt_end - peak_map.rt_start > 0.0;
        let target_rt = if peak_map.target_rt.is_finite() {
            Some(peak_map.target_rt)
        } else if has_window {
            Some((peak_map.rt_start + peak_map.rt_end) / 2.0)
        } else {
            None
        };
        let finite = |value: f64| value.is_finite().then_some(value);
        row.run_id = run_id as i64;
        row.source_file = source_file.to_owned();
        row.ms_level = ms_level;
        row.target_mz = peak_map.target_mz;
        row.target_rt = target_rt.and_then(finite);
        row.target_ion_mobility = (peak_map.target_ion_mobility >= 0.0)
            .then_some(peak_map.target_ion_mobility)
            .and_then(finite);
        row.rt_start = has_window.then_some(peak_map.rt_start).and_then(finite);
        row.rt_end = has_window.then_some(peak_map.rt_end).and_then(finite);
        row.arrays = arrays;
        if !self.pending.is_empty()
            && self.pending_bytes.saturating_add(size) > MAX_BUFFERED_BINARY_BYTES
        {
            self.flush()?;
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
        let batch = batch(&rows)?;
        let output = self
            .output
            .as_mut()
            .ok_or_else(|| Error::InvalidValue("parquet writer unavailable".into()))?;
        output.write(&batch, ROW_GROUP_ROWS)
    }

    /// Write the buffered rows and close the file (`finalize`).
    ///
    /// A second call does nothing. A consumer that received no peak map
    /// writes a file with the XIPM schema and no rows.
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

impl Drop for XIPMParquetConsumer {
    fn drop(&mut self) {
        let _ = self.finalize();
    }
}

/// The rows as a table of the XIPM schema (`buildTableFromBuilders_`).
fn batch(rows: &[Row]) -> Result<RecordBatch> {
    let mut run_id = Int64Builder::new();
    let mut source_file = StringBuilder::new();
    let mut ms_level = Int64Builder::new();
    let mut peakmap_type = StringBuilder::new();
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
    let mut target_mz = Float64Builder::new();
    let mut target_rt = Float64Builder::new();
    let mut target_im = Float64Builder::new();
    let mut rt_start = Float64Builder::new();
    let mut rt_end = Float64Builder::new();
    let mut data: [BinaryBuilder; 4] = std::array::from_fn(|_| BinaryBuilder::new());
    let mut codes: [Int64Builder; 4] = std::array::from_fn(|_| Int64Builder::new());
    let text = |builder: &mut StringBuilder, value: &str| {
        builder.append_option((!value.is_empty()).then_some(value));
    };
    for row in rows {
        run_id.append_value(row.run_id);
        text(&mut source_file, &row.source_file);
        ms_level.append_value(row.ms_level);
        peakmap_type.append_value(row.peakmap_type);
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
        target_mz.append_value(row.target_mz);
        target_rt.append_option(row.target_rt);
        target_im.append_option(row.target_ion_mobility);
        rt_start.append_option(row.rt_start);
        rt_end.append_option(row.rt_end);
        for (k, array) in row.arrays.iter().enumerate() {
            data[k].append_value(&array.bytes);
            codes[k].append_value(array.code);
        }
    }
    let [
        mut mz_data,
        mut rt_data,
        mut mobility_data,
        mut intensity_data,
    ] = data;
    let [
        mut mz_code,
        mut rt_code,
        mut mobility_code,
        mut intensity_code,
    ] = codes;
    let columns: Vec<ArrayRef> = vec![
        Arc::new(run_id.finish()),
        Arc::new(source_file.finish()),
        Arc::new(ms_level.finish()),
        Arc::new(peakmap_type.finish()),
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
        Arc::new(target_mz.finish()),
        Arc::new(target_rt.finish()),
        Arc::new(target_im.finish()),
        Arc::new(rt_start.finish()),
        Arc::new(rt_end.finish()),
        Arc::new(mz_data.finish()),
        Arc::new(rt_data.finish()),
        Arc::new(mobility_data.finish()),
        Arc::new(intensity_data.finish()),
        Arc::new(mz_code.finish()),
        Arc::new(rt_code.finish()),
        Arc::new(mobility_code.finish()),
        Arc::new(intensity_code.finish()),
    ];
    let table = RecordBatch::try_new(XIPMSchema::schema(), columns)
        .map_err(|error| Error::InvalidValue(format!("XIPM table: {error}")))?;
    let check = validation::validate(&table, &XIPMSchema::schema(), Mode::Strict);
    if !check.valid {
        return Err(Error::InvalidValue(format!(
            "XIPM table schema validation failed: {check}"
        )));
    }
    Ok(table)
}
