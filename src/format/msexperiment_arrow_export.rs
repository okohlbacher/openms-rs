// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectra and chromatograms of an [`MSExperiment`](crate::MSExperiment) as
//! Arrow tables and Parquet files.
//!
//! Ports `FORMAT/MSExperimentArrowExport.h` / `MSExperimentArrowExport.cpp`
//! (marked experimental upstream). Two layouts
//! ([`ArrowExportFormat`](crate::format::msexperiment_arrow_export::ArrowExportFormat)):
//! *long*, one row per peak with the spectrum's metadata repeated, and
//! *semi-wide*, one row per spectrum with the peaks as list columns. The
//! tables follow
//! [`SpectraLongSchema`](crate::format::arrow_schema_registry::SpectraLongSchema),
//! [`SpectraSemiWideSchema`](crate::format::arrow_schema_registry::SpectraSemiWideSchema),
//! [`ChromatogramSchema`](crate::format::arrow_schema_registry::ChromatogramSchema)
//! and
//! [`ChromatogramSemiWideSchema`](crate::format::arrow_schema_registry::ChromatogramSemiWideSchema),
//! restricted to the columns the configuration selects and the data has.
//!
//! Spectra columns (long; in semi-wide `mz`, `intensity` and `ion_mobility`
//! are lists and the scalar columns come first):
//!
//! | Column | Type | Content |
//! |---|---|---|
//! | `mz` | float64 | peak m/z |
//! | `intensity` | float32 | peak intensity |
//! | `rt` | float32 | spectrum RT (seconds) |
//! | `ion_mobility` | float32, nullable | the peak's ion mobility; only when some exported spectrum has an ion-mobility array |
//! | `spectrum_index` | uint32 | index of the spectrum in the experiment |
//! | `ms_level` | uint8 | MS level |
//! | `native_id` | utf8 | native spectrum id |
//! | `precursor_mz`, `precursor_charge` (int16), `precursor_intensity` (float32), `isolation_lower`, `isolation_upper` | nullable | the first precursor; null without one |
//!
//! Chromatogram columns: `rt` (float64), `intensity` (float32),
//! `chromatogram_index` (uint32), `native_id`, `precursor_mz`, `product_mz`.
//!
//! The source's two C Data Interface entry points hand the table to Python
//! without a copy; their counterparts here return the
//! [`RecordBatch`](arrow_array::RecordBatch) itself
//! ([`MSExperimentArrowExport::spectra_to_arrow`](crate::format::msexperiment_arrow_export::MSExperimentArrowExport::spectra_to_arrow)).
//! Exporting that batch through the C ABI needs arrow-rs's `ffi` feature,
//! which this crate does not enable yet (`docs/OPENSWATH_PARQUET_SUPPORT.md`).
//!
//! ```
//! use openms::format::msexperiment_arrow_export::{
//!     ArrowExportFormat, ArrowSpectraExportConfig, MSExperimentArrowExport,
//! };
//! use openms::{MSExperiment, MSSpectrum, Peak1D};
//!
//! let mut spectrum = MSSpectrum::default();
//! spectrum.rt = 100.0;
//! spectrum.ms_level = 1;
//! spectrum.peaks = vec![Peak1D::new(100.0, 1000.0), Peak1D::new(200.0, 2000.0)];
//! let mut experiment = MSExperiment::default();
//! experiment.spectra.push(spectrum);
//!
//! let long = MSExperimentArrowExport::spectra_to_arrow(&experiment, &ArrowSpectraExportConfig::default())?;
//! assert_eq!(long.num_rows(), 2);
//! let config = ArrowSpectraExportConfig {
//!     format: ArrowExportFormat::SemiWide,
//!     ..ArrowSpectraExportConfig::default()
//! };
//! assert_eq!(MSExperimentArrowExport::spectra_to_arrow(&experiment, &config)?.num_rows(), 1);
//! # Ok::<(), openms::Error>(())
//! ```

use crate::format::arrow_io_helpers::write_table_to_parquet;
pub use crate::format::arrow_io_helpers::{ParquetCompression, ParquetWriteConfig};
use crate::format::arrow_schema_registry::validation::{self, Mode};
use crate::format::arrow_schema_registry::{
    ChromatogramSchema, ChromatogramSemiWideSchema, SpectraLongSchema, SpectraSemiWideSchema,
};
use crate::{Error, MSExperiment, MSSpectrum, Result};
use arrow_array::builder::{
    Float32Builder, Float64Builder, Int16Builder, ListBuilder, StringBuilder, UInt8Builder,
    UInt32Builder,
};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema};
use std::collections::HashSet;
use std::path::Path;
use std::sync::Arc;

/// Layout of an Arrow export.
///
/// Ports `OpenMS::ArrowExportFormat`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ArrowExportFormat {
    /// One row per peak, the spectrum's metadata repeated (the default).
    #[default]
    Long,
    /// One row per spectrum, peaks as list columns.
    SemiWide,
}

/// Which spectra, peaks and columns a spectra export writes.
///
/// Ports `OpenMS::ArrowSpectraExportConfig`. An RT or m/z bound of exactly
/// 0 is no bound, as in the source.
#[derive(Clone, Debug, PartialEq)]
pub struct ArrowSpectraExportConfig {
    /// Layout (default [`ArrowExportFormat::Long`]).
    pub format: ArrowExportFormat,
    /// MS levels to export; empty for all.
    pub ms_levels: Vec<u32>,
    /// Minimum RT in seconds (0 for none); needs spectra sorted by RT.
    pub min_rt: f64,
    /// Maximum RT in seconds (0 for none); needs spectra sorted by RT.
    pub max_rt: f64,
    /// Minimum m/z (0 for none); needs peaks sorted by m/z.
    pub min_mz: f64,
    /// Maximum m/z (0 for none); needs peaks sorted by m/z.
    pub max_mz: f64,
    /// Columns to export by name; empty for all available.
    pub columns: Vec<String>,
    /// Export the precursor columns (default `true`).
    pub include_precursor_info: bool,
    /// Export `ion_mobility` when the data has it (default `true`).
    pub include_ion_mobility: bool,
}

impl Default for ArrowSpectraExportConfig {
    fn default() -> Self {
        Self {
            format: ArrowExportFormat::Long,
            ms_levels: Vec::new(),
            min_rt: 0.0,
            max_rt: 0.0,
            min_mz: 0.0,
            max_mz: 0.0,
            columns: Vec::new(),
            include_precursor_info: true,
            include_ion_mobility: true,
        }
    }
}

/// Which points and columns a chromatogram export writes.
///
/// Ports `OpenMS::ArrowChromatogramExportConfig`. An RT bound of exactly 0
/// is no bound.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ArrowChromatogramExportConfig {
    /// Layout (default [`ArrowExportFormat::Long`]).
    pub format: ArrowExportFormat,
    /// Minimum RT in seconds (0 for none).
    pub min_rt: f64,
    /// Maximum RT in seconds (0 for none).
    pub max_rt: f64,
    /// Columns to export by name; empty for all.
    pub columns: Vec<String>,
}

/// Export of [`MSExperiment`] data to Arrow tables and Parquet files.
///
/// Ports `OpenMS::MSExperimentArrowExport`, a class of static members.
#[derive(Clone, Copy, Debug)]
pub struct MSExperimentArrowExport;

/// `shouldIncludeColumn`.
fn wanted(name: &str, columns: &[String]) -> bool {
    columns.is_empty() || columns.iter().any(|column| column == name)
}

/// Arrow's 32-bit offsets: a utf8 column or a list column holds at most
/// `i32::MAX` bytes or elements. The source's builders fail there; arrow-rs's
/// panic, so every export counts first.
fn check_offsets(count: usize, what: &str) -> Result<()> {
    if count > i32::MAX as usize {
        return Err(Error::InvalidValue(format!(
            "{what} exceeds Arrow's 32-bit offset limit ({count} > {})",
            i32::MAX
        )));
    }
    Ok(())
}

/// The spectra the RT window selects (`getRTFilteredRange`), as an index
/// range.
fn rt_range(
    exp: &MSExperiment,
    config: &ArrowSpectraExportConfig,
) -> Result<std::ops::Range<usize>> {
    let begin = if config.min_rt != 0.0 {
        exp.rt_begin(config.min_rt)?
    } else {
        0
    };
    let end = if config.max_rt != 0.0 {
        exp.rt_end(config.max_rt)?
    } else {
        exp.spectra.len()
    };
    Ok(begin..end.max(begin))
}

/// The peaks the m/z window selects (`getMZFilteredRange`).
fn mz_range(spectrum: &MSSpectrum, min_mz: f64, max_mz: f64) -> Result<std::ops::Range<usize>> {
    let begin = if min_mz != 0.0 {
        spectrum.mz_begin(min_mz)?
    } else {
        0
    };
    let end = if max_mz != 0.0 {
        spectrum.mz_end(max_mz)?
    } else {
        spectrum.peaks.len()
    };
    Ok(begin..end.max(begin))
}

/// A selected spectrum: its index, peak range and ion-mobility values.
struct Selected<'a> {
    index: usize,
    spectrum: &'a MSSpectrum,
    peaks: std::ops::Range<usize>,
    mobility: Option<&'a [f32]>,
}

/// The spectra an export writes, in experiment order.
fn select<'a>(
    exp: &'a MSExperiment,
    config: &ArrowSpectraExportConfig,
) -> Result<Vec<Selected<'a>>> {
    let levels: HashSet<u32> = config.ms_levels.iter().copied().collect();
    let mut selected = Vec::new();
    for index in rt_range(exp, config)? {
        let spectrum = &exp.spectra[index];
        if !levels.is_empty() && !levels.contains(&spectrum.ms_level) {
            continue;
        }
        let peaks = mz_range(spectrum, config.min_mz, config.max_mz)?;
        let mobility = if spectrum.contains_im_data() {
            let (array, _) = spectrum.im_data()?;
            let values = spectrum.float_data_arrays[array].data.as_slice();
            if values.len() < peaks.end {
                return Err(Error::InvalidValue(format!(
                    "spectrum {index}: ion mobility array has {} values for {} peaks",
                    values.len(),
                    spectrum.peaks.len()
                )));
            }
            Some(values)
        } else {
            None
        };
        selected.push(Selected {
            index,
            spectrum,
            peaks,
            mobility,
        });
    }
    Ok(selected)
}

/// Whether any selected spectrum is an ion-mobility frame
/// (`experimentHasIMData`).
fn has_im_data(exp: &MSExperiment, config: &ArrowSpectraExportConfig) -> Result<bool> {
    let levels: HashSet<u32> = config.ms_levels.iter().copied().collect();
    Ok(rt_range(exp, config)?.any(|index| {
        let spectrum = &exp.spectra[index];
        (levels.is_empty() || levels.contains(&spectrum.ms_level)) && spectrum.contains_im_data()
    }))
}

/// A column: its field and finished array.
fn column(name: &str, data_type: DataType, array: ArrayRef) -> (Field, ArrayRef) {
    (Field::new(name, data_type, true), array)
}

fn table(columns: Vec<(Field, ArrayRef)>, registry: &Schema, what: &str) -> Result<RecordBatch> {
    let (fields, arrays): (Vec<Field>, Vec<ArrayRef>) = columns.into_iter().unzip();
    let schema = Arc::new(Schema::new(fields));
    let table = if arrays.is_empty() {
        RecordBatch::try_new_with_options(
            schema,
            arrays,
            &arrow_array::RecordBatchOptions::new().with_row_count(Some(0)),
        )
    } else {
        RecordBatch::try_new(schema, arrays)
    }
    .map_err(|error| Error::InvalidValue(format!("MSExperimentArrowExport: {what}: {error}")))?;
    let check = validation::validate(&table, registry, Mode::Subset);
    if !check.valid {
        return Err(Error::InvalidValue(format!(
            "MSExperimentArrowExport: {what} schema validation failed: {check}"
        )));
    }
    Ok(table)
}

/// Which precursor columns an export writes.
struct PrecursorColumns {
    mz: bool,
    charge: bool,
    intensity: bool,
    lower: bool,
    upper: bool,
}

impl PrecursorColumns {
    fn new(config: &ArrowSpectraExportConfig, names: [&str; 5]) -> Self {
        let on = |name| config.include_precursor_info && wanted(name, &config.columns);
        Self {
            mz: on(names[0]),
            charge: on(names[1]),
            intensity: on(names[2]),
            lower: on(names[3]),
            upper: on(names[4]),
        }
    }
}

/// Builders of the five precursor columns.
#[derive(Default)]
struct PrecursorBuilders {
    mz: Float64Builder,
    charge: Int16Builder,
    intensity: Float32Builder,
    lower: Float64Builder,
    upper: Float64Builder,
}

impl PrecursorBuilders {
    /// Append the first precursor of `spectrum`, or nulls.
    fn append(&mut self, spectrum: &MSSpectrum, on: &PrecursorColumns) {
        let first = spectrum.precursors.first();
        if on.mz {
            self.mz.append_option(first.map(|p| p.mz));
        }
        if on.charge {
            // `static_cast<int16_t>`: the charge wraps to 16 bits.
            self.charge.append_option(first.map(|p| p.charge as i16));
        }
        if on.intensity {
            self.intensity.append_option(first.map(|p| p.intensity));
        }
        if on.lower {
            self.lower
                .append_option(first.map(|p| p.isolation_window_lower_offset));
        }
        if on.upper {
            self.upper
                .append_option(first.map(|p| p.isolation_window_upper_offset));
        }
    }

    fn finish(mut self, on: &PrecursorColumns, names: [&str; 5], out: &mut Vec<(Field, ArrayRef)>) {
        if on.mz {
            out.push(column(
                names[0],
                DataType::Float64,
                Arc::new(self.mz.finish()),
            ));
        }
        if on.charge {
            out.push(column(
                names[1],
                DataType::Int16,
                Arc::new(self.charge.finish()),
            ));
        }
        if on.intensity {
            out.push(column(
                names[2],
                DataType::Float32,
                Arc::new(self.intensity.finish()),
            ));
        }
        if on.lower {
            out.push(column(
                names[3],
                DataType::Float64,
                Arc::new(self.lower.finish()),
            ));
        }
        if on.upper {
            out.push(column(
                names[4],
                DataType::Float64,
                Arc::new(self.upper.finish()),
            ));
        }
    }
}

const LONG_PRECURSOR: [&str; 5] = [
    SpectraLongSchema::PRECURSOR_MZ,
    SpectraLongSchema::PRECURSOR_CHARGE,
    SpectraLongSchema::PRECURSOR_INTENSITY,
    SpectraLongSchema::ISOLATION_LOWER,
    SpectraLongSchema::ISOLATION_UPPER,
];
const WIDE_PRECURSOR: [&str; 5] = [
    SpectraSemiWideSchema::PRECURSOR_MZ,
    SpectraSemiWideSchema::PRECURSOR_CHARGE,
    SpectraSemiWideSchema::PRECURSOR_INTENSITY,
    SpectraSemiWideSchema::ISOLATION_LOWER,
    SpectraSemiWideSchema::ISOLATION_UPPER,
];

/// `buildLongFormatTable`.
fn spectra_long(exp: &MSExperiment, config: &ArrowSpectraExportConfig) -> Result<RecordBatch> {
    let cols = &config.columns;
    let inc_mz = wanted(SpectraLongSchema::MZ, cols);
    let inc_intensity = wanted(SpectraLongSchema::INTENSITY, cols);
    let inc_rt = wanted(SpectraLongSchema::RT, cols);
    let inc_index = wanted(SpectraLongSchema::SPECTRUM_INDEX, cols);
    let inc_level = wanted(SpectraLongSchema::MS_LEVEL, cols);
    let inc_native = wanted(SpectraLongSchema::NATIVE_ID, cols);
    let inc_im = config.include_ion_mobility && wanted(SpectraLongSchema::ION_MOBILITY, cols);
    let precursor = PrecursorColumns::new(config, LONG_PRECURSOR);
    let has_im = inc_im && has_im_data(exp, config)?;
    let selected = select(exp, config)?;
    let mut native_bytes = 0usize;
    for s in &selected {
        native_bytes =
            native_bytes.saturating_add(s.spectrum.native_id.len().saturating_mul(s.peaks.len()));
    }
    if inc_native {
        check_offsets(native_bytes, "native_id column")?;
    }
    let mut mz = Float64Builder::new();
    let mut intensity = Float32Builder::new();
    let mut rt = Float32Builder::new();
    let mut mobility = Float32Builder::new();
    let mut index = UInt32Builder::new();
    let mut level = UInt8Builder::new();
    let mut native = StringBuilder::new();
    let mut precursors = PrecursorBuilders::default();
    for s in &selected {
        let spectrum = s.spectrum;
        for peak_index in s.peaks.clone() {
            let peak = &spectrum.peaks[peak_index];
            if inc_mz {
                mz.append_value(peak.mz);
            }
            if inc_intensity {
                intensity.append_value(peak.intensity);
            }
            if inc_rt {
                rt.append_value(spectrum.rt as f32);
            }
            if inc_index {
                // `static_cast<UInt32>`, as in the source.
                index.append_value(s.index as u32);
            }
            if inc_level {
                level.append_value(spectrum.ms_level as u8);
            }
            if inc_native {
                native.append_value(&spectrum.native_id);
            }
            if has_im {
                mobility.append_option(s.mobility.map(|values| values[peak_index]));
            }
            precursors.append(spectrum, &precursor);
        }
    }
    let mut out = Vec::new();
    if inc_mz {
        out.push(column(
            SpectraLongSchema::MZ,
            DataType::Float64,
            Arc::new(mz.finish()),
        ));
    }
    if inc_intensity {
        out.push(column(
            SpectraLongSchema::INTENSITY,
            DataType::Float32,
            Arc::new(intensity.finish()),
        ));
    }
    if inc_rt {
        out.push(column(
            SpectraLongSchema::RT,
            DataType::Float32,
            Arc::new(rt.finish()),
        ));
    }
    if has_im {
        out.push(column(
            SpectraLongSchema::ION_MOBILITY,
            DataType::Float32,
            Arc::new(mobility.finish()),
        ));
    }
    if inc_index {
        out.push(column(
            SpectraLongSchema::SPECTRUM_INDEX,
            DataType::UInt32,
            Arc::new(index.finish()),
        ));
    }
    if inc_level {
        out.push(column(
            SpectraLongSchema::MS_LEVEL,
            DataType::UInt8,
            Arc::new(level.finish()),
        ));
    }
    if inc_native {
        out.push(column(
            SpectraLongSchema::NATIVE_ID,
            DataType::Utf8,
            Arc::new(native.finish()),
        ));
    }
    precursors.finish(&precursor, LONG_PRECURSOR, &mut out);
    table(out, &SpectraLongSchema::schema(), "Spectra long format")
}

fn list_of(element: DataType) -> DataType {
    DataType::new_list(element, true)
}

/// `buildSemiWideFormatTable`.
fn spectra_semi_wide(exp: &MSExperiment, config: &ArrowSpectraExportConfig) -> Result<RecordBatch> {
    let cols = &config.columns;
    let inc_mz = wanted(SpectraSemiWideSchema::MZ, cols);
    let inc_intensity = wanted(SpectraSemiWideSchema::INTENSITY, cols);
    let inc_rt = wanted(SpectraSemiWideSchema::RT, cols);
    let inc_index = wanted(SpectraSemiWideSchema::SPECTRUM_INDEX, cols);
    let inc_level = wanted(SpectraSemiWideSchema::MS_LEVEL, cols);
    let inc_native = wanted(SpectraSemiWideSchema::NATIVE_ID, cols);
    let inc_im = config.include_ion_mobility && wanted(SpectraSemiWideSchema::ION_MOBILITY, cols);
    let precursor = PrecursorColumns::new(config, WIDE_PRECURSOR);
    let has_im = inc_im && has_im_data(exp, config)?;
    let selected = select(exp, config)?;
    let peaks: usize = selected.iter().map(|s| s.peaks.len()).sum();
    check_offsets(peaks, "peak list columns")?;
    if inc_native {
        check_offsets(
            selected.iter().map(|s| s.spectrum.native_id.len()).sum(),
            "native_id column",
        )?;
    }
    let mut index = UInt32Builder::new();
    let mut rt = Float32Builder::new();
    let mut level = UInt8Builder::new();
    let mut native = StringBuilder::new();
    let mut mz = ListBuilder::new(Float64Builder::new());
    let mut intensity = ListBuilder::new(Float32Builder::new());
    let mut mobility = ListBuilder::new(Float32Builder::new());
    let mut precursors = PrecursorBuilders::default();
    for s in &selected {
        let spectrum = s.spectrum;
        if inc_index {
            index.append_value(s.index as u32);
        }
        if inc_rt {
            rt.append_value(spectrum.rt as f32);
        }
        if inc_level {
            level.append_value(spectrum.ms_level as u8);
        }
        if inc_native {
            native.append_value(&spectrum.native_id);
        }
        precursors.append(spectrum, &precursor);
        let peaks = &spectrum.peaks[s.peaks.clone()];
        if inc_mz {
            mz.values()
                .append_slice(&peaks.iter().map(|p| p.mz).collect::<Vec<_>>());
            mz.append(true);
        }
        if inc_intensity {
            intensity
                .values()
                .append_slice(&peaks.iter().map(|p| p.intensity).collect::<Vec<_>>());
            intensity.append(true);
        }
        if has_im {
            match s.mobility {
                Some(values) => {
                    mobility.values().append_slice(&values[s.peaks.clone()]);
                    mobility.append(true);
                }
                None => mobility.append(false),
            }
        }
    }
    let mut out = Vec::new();
    if inc_index {
        out.push(column(
            SpectraSemiWideSchema::SPECTRUM_INDEX,
            DataType::UInt32,
            Arc::new(index.finish()),
        ));
    }
    if inc_rt {
        out.push(column(
            SpectraSemiWideSchema::RT,
            DataType::Float32,
            Arc::new(rt.finish()),
        ));
    }
    if inc_level {
        out.push(column(
            SpectraSemiWideSchema::MS_LEVEL,
            DataType::UInt8,
            Arc::new(level.finish()),
        ));
    }
    if inc_native {
        out.push(column(
            SpectraSemiWideSchema::NATIVE_ID,
            DataType::Utf8,
            Arc::new(native.finish()),
        ));
    }
    if inc_mz {
        out.push(column(
            SpectraSemiWideSchema::MZ,
            list_of(DataType::Float64),
            Arc::new(mz.finish()),
        ));
    }
    if inc_intensity {
        out.push(column(
            SpectraSemiWideSchema::INTENSITY,
            list_of(DataType::Float32),
            Arc::new(intensity.finish()),
        ));
    }
    if has_im {
        out.push(column(
            SpectraSemiWideSchema::ION_MOBILITY,
            list_of(DataType::Float32),
            Arc::new(mobility.finish()),
        ));
    }
    precursors.finish(&precursor, WIDE_PRECURSOR, &mut out);
    table(
        out,
        &SpectraSemiWideSchema::schema(),
        "Spectra semi-wide format",
    )
}

/// Whether a chromatogram point passes the RT window.
fn rt_passes(rt: f64, config: &ArrowChromatogramExportConfig) -> bool {
    !(config.min_rt != 0.0 && rt < config.min_rt) && !(config.max_rt != 0.0 && rt > config.max_rt)
}

/// `exportChromatogramsToArrow`.
fn chromatograms(
    exp: &MSExperiment,
    config: &ArrowChromatogramExportConfig,
) -> Result<RecordBatch> {
    let cols = &config.columns;
    let inc_rt = wanted(ChromatogramSchema::RT, cols);
    let inc_intensity = wanted(ChromatogramSchema::INTENSITY, cols);
    let inc_index = wanted(ChromatogramSchema::CHROMATOGRAM_INDEX, cols);
    let inc_native = wanted(ChromatogramSchema::NATIVE_ID, cols);
    let inc_precursor = wanted(ChromatogramSchema::PRECURSOR_MZ, cols);
    let inc_product = wanted(ChromatogramSchema::PRODUCT_MZ, cols);
    let mut index = UInt32Builder::new();
    let mut native = StringBuilder::new();
    let mut precursor = Float64Builder::new();
    let mut product = Float64Builder::new();
    let mut out = Vec::new();
    if config.format == ArrowExportFormat::Long {
        let mut native_bytes = 0usize;
        for chromatogram in &exp.chromatograms {
            let kept = chromatogram
                .peaks
                .iter()
                .filter(|point| rt_passes(point.rt, config))
                .count();
            native_bytes =
                native_bytes.saturating_add(chromatogram.native_id.len().saturating_mul(kept));
        }
        if inc_native {
            check_offsets(native_bytes, "native_id column")?;
        }
        let mut rt = Float64Builder::new();
        let mut intensity = Float32Builder::new();
        for (chrom_index, chromatogram) in exp.chromatograms.iter().enumerate() {
            for point in chromatogram
                .peaks
                .iter()
                .filter(|point| rt_passes(point.rt, config))
            {
                if inc_rt {
                    rt.append_value(point.rt);
                }
                if inc_intensity {
                    intensity.append_value(point.intensity);
                }
                if inc_index {
                    index.append_value(chrom_index as u32);
                }
                if inc_native {
                    native.append_value(&chromatogram.native_id);
                }
                if inc_precursor {
                    precursor.append_value(chromatogram.precursor.mz);
                }
                if inc_product {
                    product.append_value(chromatogram.product.mz);
                }
            }
        }
        if inc_rt {
            out.push(column(
                ChromatogramSchema::RT,
                DataType::Float64,
                Arc::new(rt.finish()),
            ));
        }
        if inc_intensity {
            out.push(column(
                ChromatogramSchema::INTENSITY,
                DataType::Float32,
                Arc::new(intensity.finish()),
            ));
        }
        if inc_index {
            out.push(column(
                ChromatogramSchema::CHROMATOGRAM_INDEX,
                DataType::UInt32,
                Arc::new(index.finish()),
            ));
        }
        if inc_native {
            out.push(column(
                ChromatogramSchema::NATIVE_ID,
                DataType::Utf8,
                Arc::new(native.finish()),
            ));
        }
        if inc_precursor {
            out.push(column(
                ChromatogramSchema::PRECURSOR_MZ,
                DataType::Float64,
                Arc::new(precursor.finish()),
            ));
        }
        if inc_product {
            out.push(column(
                ChromatogramSchema::PRODUCT_MZ,
                DataType::Float64,
                Arc::new(product.finish()),
            ));
        }
        return table(
            out,
            &ChromatogramSchema::schema(),
            "Chromatogram long format",
        );
    }
    let points: usize = exp.chromatograms.iter().map(|c| c.peaks.len()).sum();
    check_offsets(points, "point list columns")?;
    if inc_native {
        check_offsets(
            exp.chromatograms.iter().map(|c| c.native_id.len()).sum(),
            "native_id column",
        )?;
    }
    let mut rt = ListBuilder::new(Float64Builder::new());
    let mut intensity = ListBuilder::new(Float32Builder::new());
    for (chrom_index, chromatogram) in exp.chromatograms.iter().enumerate() {
        if inc_index {
            index.append_value(chrom_index as u32);
        }
        if inc_native {
            native.append_value(&chromatogram.native_id);
        }
        if inc_precursor {
            precursor.append_value(chromatogram.precursor.mz);
        }
        if inc_product {
            product.append_value(chromatogram.product.mz);
        }
        let kept = chromatogram
            .peaks
            .iter()
            .filter(|point| rt_passes(point.rt, config));
        if inc_rt {
            for point in kept.clone() {
                rt.values().append_value(point.rt);
            }
            rt.append(true);
        }
        if inc_intensity {
            for point in kept {
                intensity.values().append_value(point.intensity);
            }
            intensity.append(true);
        }
    }
    let semi = ChromatogramSemiWideSchema::CHROMATOGRAM_INDEX;
    if inc_index {
        out.push(column(semi, DataType::UInt32, Arc::new(index.finish())));
    }
    if inc_native {
        out.push(column(
            ChromatogramSemiWideSchema::NATIVE_ID,
            DataType::Utf8,
            Arc::new(native.finish()),
        ));
    }
    if inc_rt {
        out.push(column(
            ChromatogramSemiWideSchema::RT,
            list_of(DataType::Float64),
            Arc::new(rt.finish()),
        ));
    }
    if inc_intensity {
        out.push(column(
            ChromatogramSemiWideSchema::INTENSITY,
            list_of(DataType::Float32),
            Arc::new(intensity.finish()),
        ));
    }
    if inc_precursor {
        out.push(column(
            ChromatogramSemiWideSchema::PRECURSOR_MZ,
            DataType::Float64,
            Arc::new(precursor.finish()),
        ));
    }
    if inc_product {
        out.push(column(
            ChromatogramSemiWideSchema::PRODUCT_MZ,
            DataType::Float64,
            Arc::new(product.finish()),
        ));
    }
    table(
        out,
        &ChromatogramSemiWideSchema::schema(),
        "Chromatogram semi-wide format",
    )
}

fn filter_names(names: Vec<&str>, columns: &[String]) -> Vec<String> {
    names
        .into_iter()
        .filter(|name| wanted(name, columns))
        .map(str::to_owned)
        .collect()
}

impl MSExperimentArrowExport {
    /// The columns a spectra export would write (`getSpectraArrowColumnNames`),
    /// in table order: the layout's columns, `ion_mobility` only when some
    /// selected spectrum has ion-mobility data, the precursor columns only
    /// with [`include_precursor_info`](ArrowSpectraExportConfig::include_precursor_info),
    /// then restricted to [`columns`](ArrowSpectraExportConfig::columns) when
    /// that is not empty.
    ///
    /// # Errors
    ///
    /// [`Error::UnsortedData`] when an RT bound is set and the spectra are not
    /// sorted by RT (the source's binary search assumes it silently).
    pub fn spectra_arrow_column_names(
        exp: &MSExperiment,
        config: &ArrowSpectraExportConfig,
    ) -> Result<Vec<String>> {
        let has_im = config.include_ion_mobility && has_im_data(exp, config)?;
        let mut names: Vec<&str> = match config.format {
            ArrowExportFormat::Long => {
                let mut names = vec![
                    SpectraLongSchema::MZ,
                    SpectraLongSchema::INTENSITY,
                    SpectraLongSchema::RT,
                ];
                if has_im {
                    names.push(SpectraLongSchema::ION_MOBILITY);
                }
                names.extend([
                    SpectraLongSchema::SPECTRUM_INDEX,
                    SpectraLongSchema::MS_LEVEL,
                    SpectraLongSchema::NATIVE_ID,
                ]);
                names
            }
            ArrowExportFormat::SemiWide => {
                let mut names = vec![
                    SpectraSemiWideSchema::SPECTRUM_INDEX,
                    SpectraSemiWideSchema::RT,
                    SpectraSemiWideSchema::MS_LEVEL,
                    SpectraSemiWideSchema::NATIVE_ID,
                    SpectraSemiWideSchema::MZ,
                    SpectraSemiWideSchema::INTENSITY,
                ];
                if has_im {
                    names.push(SpectraSemiWideSchema::ION_MOBILITY);
                }
                names
            }
        };
        if config.include_precursor_info {
            names.extend(match config.format {
                ArrowExportFormat::Long => LONG_PRECURSOR,
                ArrowExportFormat::SemiWide => WIDE_PRECURSOR,
            });
        }
        Ok(filter_names(names, &config.columns))
    }

    /// The columns a chromatogram export would write
    /// (`getChromatogramArrowColumnNames`); the experiment is not consulted,
    /// as in the source.
    pub fn chromatogram_arrow_column_names(
        _exp: &MSExperiment,
        config: &ArrowChromatogramExportConfig,
    ) -> Vec<String> {
        let names = match config.format {
            ArrowExportFormat::Long => vec![
                ChromatogramSchema::RT,
                ChromatogramSchema::INTENSITY,
                ChromatogramSchema::CHROMATOGRAM_INDEX,
                ChromatogramSchema::NATIVE_ID,
                ChromatogramSchema::PRECURSOR_MZ,
                ChromatogramSchema::PRODUCT_MZ,
            ],
            ArrowExportFormat::SemiWide => vec![
                ChromatogramSemiWideSchema::CHROMATOGRAM_INDEX,
                ChromatogramSemiWideSchema::NATIVE_ID,
                ChromatogramSemiWideSchema::RT,
                ChromatogramSemiWideSchema::INTENSITY,
                ChromatogramSemiWideSchema::PRECURSOR_MZ,
                ChromatogramSemiWideSchema::PRODUCT_MZ,
            ],
        };
        filter_names(names, &config.columns)
    }

    /// The spectra as one Arrow table (the table the source's
    /// `exportSpectraToArrowCDataInterface` hands over).
    ///
    /// Spectra are selected by the RT window (a binary search, so the
    /// spectra must be sorted by RT when a bound is set) and the MS levels,
    /// peaks by the m/z window (sorted peaks needed when a bound is set).
    /// `spectrum_index` is the spectrum's index in the whole experiment. An
    /// experiment with nothing selected gives a table with the columns and no
    /// rows.
    ///
    /// # Errors
    ///
    /// [`Error::UnsortedData`] as above; [`Error::InvalidValue`] for a
    /// non-finite bound, an ion-mobility array shorter than its spectrum's
    /// peak list (the source reads past its end), or a column beyond Arrow's
    /// 32-bit offsets (2^31 - 1 string bytes or list elements; the source's
    /// builder fails there).
    pub fn spectra_to_arrow(
        exp: &MSExperiment,
        config: &ArrowSpectraExportConfig,
    ) -> Result<RecordBatch> {
        match config.format {
            ArrowExportFormat::Long => spectra_long(exp, config),
            ArrowExportFormat::SemiWide => spectra_semi_wide(exp, config),
        }
    }

    /// The chromatograms as one Arrow table (the table the source's
    /// `exportChromatogramsToArrowCDataInterface` hands over).
    ///
    /// Points are selected by the RT window point by point (no sorting
    /// needed). `chromatogram_index` counts every chromatogram, whether or
    /// not any of its points is kept; in the semi-wide layout a chromatogram
    /// without kept points is a row with empty lists.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a column beyond Arrow's 32-bit offsets.
    pub fn chromatograms_to_arrow(
        exp: &MSExperiment,
        config: &ArrowChromatogramExportConfig,
    ) -> Result<RecordBatch> {
        chromatograms(exp, config)
    }

    /// Write the spectra table to a Parquet file (`exportSpectraToParquet`),
    /// with [`write_table_to_parquet`]: the Arrow schema is stored, and
    /// `parquet_config` selects codec, level, row group size, statistics and
    /// page size. The destination is replaced only when the whole file has
    /// been written, so a failed export never leaves a partial file (the
    /// source removes its fragment instead).
    ///
    /// # Errors
    ///
    /// As [`spectra_to_arrow`](Self::spectra_to_arrow) and
    /// [`write_table_to_parquet`] (for example a row group size of 0 for a
    /// table with rows).
    pub fn export_spectra_to_parquet(
        exp: &MSExperiment,
        filename: impl AsRef<Path>,
        config: &ArrowSpectraExportConfig,
        parquet_config: &ParquetWriteConfig,
    ) -> Result<()> {
        let table = Self::spectra_to_arrow(exp, config)?;
        write_table_to_parquet(&table, filename, parquet_config)
    }

    /// Write the chromatogram table to a Parquet file
    /// (`exportChromatogramsToParquet`); see
    /// [`export_spectra_to_parquet`](Self::export_spectra_to_parquet).
    ///
    /// # Errors
    ///
    /// As [`chromatograms_to_arrow`](Self::chromatograms_to_arrow) and
    /// [`write_table_to_parquet`].
    pub fn export_chromatograms_to_parquet(
        exp: &MSExperiment,
        filename: impl AsRef<Path>,
        config: &ArrowChromatogramExportConfig,
        parquet_config: &ParquetWriteConfig,
    ) -> Result<()> {
        let table = Self::chromatograms_to_arrow(exp, config)?;
        write_table_to_parquet(&table, filename, parquet_config)
    }
}
