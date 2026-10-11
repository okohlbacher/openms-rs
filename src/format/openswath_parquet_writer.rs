// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Shared machinery of the OpenSWATH Parquet writers.
//!
//! The source's `MSChromatogramParquetConsumer.cpp`,
//! `MobilogramParquetConsumer.cpp` and `XIPMParquetConsumer.cpp` each repeat
//! the same three pieces in an anonymous namespace: the transition-library
//! lookup tables (`buildTransitionMaps_`, `parseOrAssignId_`,
//! `buildAnnotation_`), the Numpress + zlib array encoding, and a Parquet
//! `FileWriter` with ZSTD level 11 that is fed chunk by chunk. This module
//! holds one copy of each, parameterised where the three files differ.

use crate::analysis::openswath::transition_experiment::LightTargetedExperiment;
use crate::format::numpress_coder::{
    self as codec, NumpressCoderLimits, NumpressCompression, NumpressConfig, NumpressEncodeStatus,
    Work,
};
use crate::format::parquet_reader_helpers::source_to_int64;
use crate::format::path_io::StagedFile;
use crate::{Error, Result};
use arrow_array::RecordBatch;
use arrow_schema::SchemaRef;
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_writer::ArrowWriterOptions;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};

/// Precursor-level metadata of a compound (`CompoundInfo`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct CompoundInfo {
    pub(crate) precursor_id: i64,
    pub(crate) modified_sequence: String,
    pub(crate) precursor_charge: i64,
    pub(crate) precursor_decoy: i64,
}

/// Transition-level metadata (`TransitionInfo`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct TransitionInfo {
    pub(crate) transition_id: i64,
    pub(crate) precursor_id: i64,
    pub(crate) modified_sequence: String,
    pub(crate) precursor_charge: i64,
    pub(crate) product_charge: i64,
    pub(crate) detecting_transition: i64,
    pub(crate) precursor_decoy: i64,
    pub(crate) product_decoy: i64,
    pub(crate) transition_ordinal: i64,
    pub(crate) transition_type: String,
    pub(crate) annotation: String,
}

/// The lookup tables `buildTransitionMaps_` builds from a transition library.
#[derive(Clone, Debug, Default)]
pub(crate) struct TransitionMaps {
    /// By compound id (or, for a transition whose `peptide_ref` names no
    /// compound, by that reference).
    pub(crate) compounds: HashMap<String, CompoundInfo>,
    /// By transition name.
    pub(crate) transitions: HashMap<String, TransitionInfo>,
}

/// `parseOrAssignId_`: the id a library text maps to.
///
/// A text that `StringUtils::toInt64` accepts, and that is not yet taken, is
/// its own id (and moves the counter past it); anything else gets the lowest
/// free id at or above the counter. `reject_negative` is the mobilogram
/// writer's extra rule that a negative parsed id is never kept.
fn parse_or_assign_id(
    text: &str,
    next_id: &mut i64,
    used: &mut HashSet<i64>,
    reject_negative: bool,
) -> i64 {
    if let Some(value) = source_to_int64(text) {
        if (!reject_negative || value >= 0) && used.insert(value) {
            if value >= *next_id {
                *next_id = value.saturating_add(1);
            }
            return value;
        }
    }
    while used.contains(next_id) {
        *next_id = next_id.saturating_add(1);
    }
    used.insert(*next_id);
    let id = *next_id;
    *next_id = next_id.saturating_add(1);
    id
}

/// `buildAnnotation_`: `<type><ordinal>` plus `^<charge>` for a positive
/// charge; empty without a type or with a negative ordinal.
fn build_annotation(transition_type: &str, ordinal: i64, charge: i64) -> String {
    if transition_type.is_empty() || ordinal < 0 {
        return String::new();
    }
    let mut annotation = format!("{transition_type}{ordinal}");
    if charge > 0 {
        annotation.push_str(&format!("^{charge}"));
    }
    annotation
}

impl TransitionMaps {
    /// `buildTransitionMaps_`.
    ///
    /// Compounds are numbered first, in library order, then transitions. A
    /// precursor is a decoy only when it has decoy transitions and no target
    /// transition (IPF libraries mix target detecting and decoy identifying
    /// transitions of one precursor). A duplicated compound id or
    /// transition name keeps its first entry, as the source's `emplace` does;
    /// a duplicated compound id still consumes a second id, which the
    /// transitions then reference, as in the source.
    pub(crate) fn build(library: &LightTargetedExperiment, reject_negative_ids: bool) -> Self {
        let mut has_target: HashSet<&str> = HashSet::new();
        let mut has_decoy: HashSet<&str> = HashSet::new();
        for transition in &library.transitions {
            if transition.decoy() {
                has_decoy.insert(&transition.peptide_ref);
            } else {
                has_target.insert(&transition.peptide_ref);
            }
        }
        let decoy_of = |reference: &str| {
            i64::from(!has_target.contains(reference) && has_decoy.contains(reference))
        };
        let mut next_precursor_id = 1;
        let mut next_transition_id = 1;
        let mut used_precursor_ids = HashSet::new();
        let mut used_transition_ids = HashSet::new();
        let mut maps = Self::default();
        let mut precursor_ids: HashMap<String, i64> = HashMap::new();
        for compound in &library.compounds {
            let precursor_id = parse_or_assign_id(
                &compound.id,
                &mut next_precursor_id,
                &mut used_precursor_ids,
                reject_negative_ids,
            );
            precursor_ids.insert(compound.id.clone(), precursor_id);
            maps.compounds
                .entry(compound.id.clone())
                .or_insert_with(|| CompoundInfo {
                    precursor_id,
                    modified_sequence: compound.sequence.clone(),
                    precursor_charge: i64::from(compound.charge),
                    precursor_decoy: decoy_of(&compound.id),
                });
        }
        let mut transition_ids: HashMap<String, i64> = HashMap::new();
        for transition in &library.transitions {
            let name = &transition.transition_name;
            let reference = &transition.peptide_ref;
            let transition_id = match transition_ids.get(name) {
                Some(&id) => id,
                None => {
                    let id = parse_or_assign_id(
                        name,
                        &mut next_transition_id,
                        &mut used_transition_ids,
                        reject_negative_ids,
                    );
                    transition_ids.insert(name.clone(), id);
                    id
                }
            };
            let precursor_id = match precursor_ids.get(reference) {
                Some(&id) => id,
                None => {
                    let id = parse_or_assign_id(
                        reference,
                        &mut next_precursor_id,
                        &mut used_precursor_ids,
                        reject_negative_ids,
                    );
                    precursor_ids.insert(reference.clone(), id);
                    maps.compounds
                        .entry(reference.clone())
                        .or_insert_with(|| CompoundInfo {
                            precursor_id: id,
                            modified_sequence: String::new(),
                            precursor_charge: 0,
                            precursor_decoy: decoy_of(reference),
                        });
                    id
                }
            };
            let mut info = TransitionInfo {
                transition_id,
                precursor_id,
                ..TransitionInfo::default()
            };
            if let Some(compound) = maps.compounds.get(reference) {
                info.modified_sequence = compound.modified_sequence.clone();
                info.precursor_charge = compound.precursor_charge;
                info.precursor_decoy = compound.precursor_decoy;
            }
            info.product_charge = i64::from(transition.fragment_charge);
            info.detecting_transition = i64::from(transition.is_detecting_transition());
            info.product_decoy = i64::from(transition.decoy());
            info.transition_ordinal = i64::from(transition.fragment_nr);
            info.transition_type = transition.fragment_type_str().to_owned();
            info.annotation = build_annotation(
                &info.transition_type,
                info.transition_ordinal,
                info.product_charge,
            );
            maps.transitions.entry(name.clone()).or_insert(info);
        }
        maps
    }
}

/// Compression identifier 1: zlib-compressed little-endian doubles.
pub(crate) const ZLIB_DOUBLES: i64 = 1;
/// Compression identifier 5: MS-Numpress linear, then zlib.
pub(crate) const NUMPRESS_LINEAR_ZLIB: i64 = 5;
/// Compression identifier 6: MS-Numpress short logged float, then zlib.
pub(crate) const NUMPRESS_SLOF_ZLIB: i64 = 6;

/// How one array is encoded lossily, as the source configures its
/// `NumpressConfig`.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Lossy {
    /// MS-Numpress linear with an estimated fixed point; `mass_accuracy` is
    /// `linear_fp_mass_acc` (the chromatogram writer's 0.05, or `None` for
    /// the maximal fixed point).
    Linear { mass_accuracy: Option<f64> },
    /// MS-Numpress short logged float with an estimated fixed point.
    Slof,
}

/// One encoded array: its compression identifier and bytes.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Encoded {
    pub(crate) code: i64,
    pub(crate) bytes: Vec<u8>,
}

fn zlib(bytes: &[u8], limits: NumpressCoderLimits) -> Result<Vec<u8>> {
    codec::zlib_encode(bytes, &mut Work::new(limits))
}

/// The lossless fallback: identifier 1, zlib over little-endian doubles.
fn lossless(values: &[f64], limits: NumpressCoderLimits) -> Result<Encoded> {
    let mut raw = Vec::with_capacity(values.len().saturating_mul(8));
    for value in values {
        raw.extend_from_slice(&value.to_le_bytes());
    }
    Ok(Encoded {
        code: ZLIB_DOUBLES,
        bytes: zlib(&raw, limits)?,
    })
}

/// Encode one array the way the source does — MS-Numpress, then zlib —
/// unless that would not round-trip.
///
/// An empty array is empty bytes under the lossy identifier, as in the
/// source. Native difference: where the source's Numpress output cannot be
/// decoded back to the data — a linear fixed point that is not positive
/// (the mass-accuracy estimator returns 0 for fewer than three values, and
/// the source then stores every value as 0, which decodes as NaN), first
/// values that do not quantise into the unsigned 32-bit range linear stores
/// them in (negative or huge retention times wrap), negative or non-finite
/// intensities for slof, or an encoder rejection — the array is stored
/// losslessly under identifier 1, which every source reader decodes.
pub(crate) fn encode(values: &[f64], lossy: Lossy, limits: NumpressCoderLimits) -> Result<Encoded> {
    let code = match lossy {
        Lossy::Linear { .. } => NUMPRESS_LINEAR_ZLIB,
        Lossy::Slof => NUMPRESS_SLOF_ZLIB,
    };
    if values.is_empty() {
        return Ok(Encoded {
            code,
            bytes: Vec::new(),
        });
    }
    let config = match lossy {
        Lossy::Linear { mass_accuracy } => NumpressConfig {
            compression: NumpressCompression::Linear,
            estimate_fixed_point: true,
            error_tolerance: -1.0,
            linear_fp_mass_acc: mass_accuracy.unwrap_or(-1.0),
            ..NumpressConfig::default()
        },
        Lossy::Slof => NumpressConfig {
            compression: NumpressCompression::Slof,
            estimate_fixed_point: true,
            error_tolerance: -1.0,
            ..NumpressConfig::default()
        },
    };
    let usable_input = match lossy {
        Lossy::Linear { .. } => values.iter().all(|value| value.is_finite()),
        Lossy::Slof => values
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0),
    };
    if !usable_input {
        return lossless(values, limits);
    }
    let report = codec::encode_raw(values, &config, &mut Work::new(limits))?;
    if report.status != NumpressEncodeStatus::Encoded {
        return lossless(values, limits);
    }
    let factor = report.fixed_point.unwrap_or(0.0);
    if !(factor.is_finite() && factor > 0.0) {
        return lossless(values, limits);
    }
    if matches!(lossy, Lossy::Linear { .. })
        && values.iter().take(2).any(|value| {
            let quantized = (value * factor + 0.5).trunc();
            !(0.0..4_294_967_296.0).contains(&quantized)
        })
    {
        return lossless(values, limits);
    }
    Ok(Encoded {
        code,
        bytes: zlib(&report.output, limits)?,
    })
}

/// A Parquet file written row group by row group through a staged
/// temporary file: the source's `parquet::arrow::FileWriter` with ZSTD level
/// 11 and no stored Arrow schema.
///
/// Native difference: the source opens (and truncates) the destination at
/// its first flush and leaves a partial file when a later write fails; this
/// writes a sibling temporary file and moves it onto the destination only in
/// [`finish`](Self::finish).
pub(crate) struct StagedParquet {
    destination: PathBuf,
    schema: SchemaRef,
    open: Option<(StagedFile, ArrowWriter<File>)>,
    rows_written: usize,
}

impl StagedParquet {
    pub(crate) fn new(destination: &Path, schema: SchemaRef) -> Self {
        Self {
            destination: destination.to_path_buf(),
            schema,
            open: None,
            rows_written: 0,
        }
    }

    fn writer(&mut self) -> Result<&mut ArrowWriter<File>> {
        if self.open.is_none() {
            let (staged, file) = StagedFile::create(&self.destination)?;
            let level = ZstdLevel::try_new(11)
                .map_err(|error| Error::InvalidValue(format!("ZSTD level: {error}")))?;
            let properties = WriterProperties::builder()
                .set_compression(Compression::ZSTD(level))
                .build();
            let options = ArrowWriterOptions::new()
                .with_properties(properties)
                .with_skip_arrow_metadata(true);
            let writer = ArrowWriter::try_new_with_options(file, self.schema.clone(), options)
                .map_err(|error| {
                    Error::InvalidValue(format!("Failed to create parquet FileWriter: {error}"))
                })?;
            self.open = Some((staged, writer));
        }
        match self.open.as_mut() {
            Some((_, writer)) => Ok(writer),
            None => Err(Error::InvalidValue("parquet writer unavailable".into())),
        }
    }

    /// Write `batch` as row groups of at most `row_group_rows` rows, as the
    /// source's `WriteTable(table, chunk_size)` does.
    pub(crate) fn write(&mut self, batch: &RecordBatch, row_group_rows: usize) -> Result<()> {
        let row_group_rows = row_group_rows.max(1);
        let writer = self.writer()?;
        let fail = |error: parquet::errors::ParquetError| {
            Error::InvalidValue(format!("Failed to write parquet table: {error}"))
        };
        let mut offset = 0;
        while offset < batch.num_rows() {
            let length = row_group_rows.min(batch.num_rows() - offset);
            writer.write(&batch.slice(offset, length)).map_err(fail)?;
            writer.flush().map_err(fail)?;
            offset += length;
        }
        self.rows_written += batch.num_rows();
        Ok(())
    }

    /// Close the file and move it onto the destination. A file that received
    /// no row gets one empty row group, as the source's `WriteTable` of an
    /// empty table writes one.
    pub(crate) fn finish(mut self) -> Result<()> {
        let fail = |error: parquet::errors::ParquetError| {
            Error::InvalidValue(format!("Failed to close parquet writer: {error}"))
        };
        if self.rows_written == 0 {
            let writer = match self.open.take() {
                Some((staged, writer)) => (staged, writer),
                None => {
                    self.writer()?;
                    self.open
                        .take()
                        .ok_or_else(|| Error::InvalidValue("parquet writer unavailable".into()))?
                }
            };
            let (staged, writer) = writer;
            let (mut file_writer, factory) = writer.into_serialized_writer().map_err(fail)?;
            let columns = factory.create_column_writers(0).map_err(fail)?;
            let mut group = file_writer.next_row_group().map_err(fail)?;
            for column in columns {
                column
                    .close()
                    .map_err(fail)?
                    .append_to_row_group(&mut group)
                    .map_err(fail)?;
            }
            group.close().map_err(fail)?;
            let file = file_writer.into_inner().map_err(fail)?;
            return staged.commit(file, &self.destination);
        }
        let (staged, writer) = self
            .open
            .take()
            .ok_or_else(|| Error::InvalidValue("parquet writer unavailable".into()))?;
        let file = writer.into_inner().map_err(fail)?;
        staged.commit(file, &self.destination)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_follow_parse_or_assign() {
        let mut next = 1;
        let mut used = HashSet::new();
        assert_eq!(parse_or_assign_id("5", &mut next, &mut used, false), 5);
        assert_eq!(next, 6);
        assert_eq!(parse_or_assign_id("pep", &mut next, &mut used, false), 6);
        assert_eq!(parse_or_assign_id("5", &mut next, &mut used, false), 7);
        assert_eq!(parse_or_assign_id(" 2 ", &mut next, &mut used, false), 2);
        assert_eq!(parse_or_assign_id("-3", &mut next, &mut used, false), -3);
        assert_eq!(parse_or_assign_id("-4", &mut next, &mut used, true), 8);
        assert_eq!(build_annotation("y", 7, 1), "y7^1");
        assert_eq!(build_annotation("b", 3, 0), "b3");
        assert_eq!(build_annotation("", 3, 1), "");
        assert_eq!(build_annotation("y", -1, 1), "");
    }
}
