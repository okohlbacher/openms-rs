// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned spectra, chromatograms and experiments.
//!
//! Peak coordinates use `f64`; intensities use `f32`, matching OpenMS. Fields are
//! public for ergonomic construction. Call [`MSSpectrum::validate`] after direct
//! edits; checked sorting and selection preserve all parallel data arrays.
//! Searches validate coordinate order on each call (O(n)), then perform a binary
//! search. Empty searches return `None`, and ranges are recomputed on demand so
//! that public mutation cannot leave a stale range cache.

use crate::error::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::hash::{Hash, Hasher};

mod area_iteration;
mod experiment_2d;
mod peak2d;
pub use experiment_2d::Data2DLimits;
pub use peak2d::{MobilityPeak2D, Peak2D, RichPeak2D};
mod peak_data;
mod peak_index;
pub use area_iteration::{
    AreaBounds, AreaIter, AreaIterMut, AreaLimits, AreaOptions, AreaPeak, AreaPeakMut,
};
pub use peak_data::{FlatPeakData, PeakDataLimits, SpectrumPeakData};
pub use peak_index::PeakIndex;
mod acquisition_fields;
mod chromatogram_tools;
mod data_array;
mod mass_trace;
pub use chromatogram_tools::{
    ChromatogramConversionLimits, ChromatogramConversionReport, ChromatogramTools,
};
pub use mass_trace::{MassTrace, MassTraceLimits, MassTraceQuantMethod};
mod experiment_aggregation;
mod experiment_summary;
mod mobilogram;
pub use experiment_summary::SummaryLimits;
pub use mobilogram::{MobilityPeak1D, Mobilogram, MobilogramLimits, MobilogramRanges};
pub mod features;
pub mod geometry;
pub use experiment_aggregation::{AggregationLimits, MzAggregation, MzRtRegion};
pub use features::{
    BaseFeature, ColumnHeader, ConsensusFeature, ConsensusMap, Feature, FeatureHandle, FeatureMap,
    FeatureRanges,
};
pub use geometry::{BoundingBox2D, ConvexHull2D, Point2D};

/// A single mass-to-charge measurement and its intensity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Peak1D {
    pub mz: f64,
    pub intensity: f32,
}

impl Peak1D {
    pub const fn new(mz: f64, intensity: f32) -> Self {
        Self { mz, intensity }
    }
}

/// A chromatographic measurement, with retention time in seconds.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChromatogramPeak {
    pub rt: f64,
    pub intensity: f32,
}

impl ChromatogramPeak {
    pub const fn new(rt: f64, intensity: f32) -> Self {
        Self { rt, intensity }
    }
}

// The source's comparator overloads map directly to scalar comparisons on the
// public fields. Formatting and hashing are the remaining value operations.
macro_rules! peak_value_traits {
    ($peak:ty, $position:ident) => {
        impl fmt::Display for $peak {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                match formatter.precision() {
                    Some(precision) => write!(
                        formatter,
                        "POS: {:.*} INT: {:.*}",
                        precision, self.$position, precision, self.intensity
                    ),
                    None => write!(formatter, "POS: {} INT: {}", self.$position, self.intensity),
                }
            }
        }

        impl Hash for $peak {
            fn hash<H: Hasher>(&self, state: &mut H) {
                // Ordinary floating equality identifies the two zero signs.
                // Other bits, including NaN payloads, remain distinct inputs.
                let position = if self.$position == 0.0 {
                    0
                } else {
                    self.$position.to_bits()
                };
                let intensity = if self.intensity == 0.0 {
                    0
                } else {
                    self.intensity.to_bits()
                };
                position.hash(state);
                intensity.hash(state);
            }
        }
    };
}

peak_value_traits!(Peak1D, mz);
peak_value_traits!(ChromatogramPeak, rt);

/// Selected precursor ion and its acquisition metadata. Charge zero means unknown.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Precursor {
    pub mz: f64,
    pub intensity: f32,
    pub charge: i32,
    pub activation_methods: BTreeSet<crate::metadata::ActivationMethod>,
    pub activation_energy: f64,
    pub isolation_window_lower_offset: f64,
    pub isolation_window_upper_offset: f64,
    /// An isolation target different from the selected ion m/z; None uses mz.
    pub isolation_target_mz: Option<f64>,
    pub drift_time: Option<f64>,
    pub drift_time_unit: crate::metadata::DriftTimeUnit,
    pub drift_window_lower_offset: f64,
    pub drift_window_upper_offset: f64,
    pub possible_charge_states: Vec<i32>,
    pub cv_terms: crate::metadata::CVTermList,
    /// Native ID of the parent spectrum, corresponding to mzML spectrumRef.
    pub spectrum_reference: Option<String>,
}

impl Precursor {
    pub fn new(mz: f64, charge: i32) -> Self {
        Self {
            mz,
            intensity: 0.0,
            charge,
            ..Self::default()
        }
    }

    pub(crate) fn has_acquisition_metadata(&self) -> bool {
        let mut basic = Self::new(self.mz, self.charge);
        basic.intensity = self.intensity;
        self != &basic
    }
}

/// Named per-peak values. An empty array is permitted as a placeholder.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DataArray<T> {
    pub name: String,
    pub data: Vec<T>,
    /// Source MetaInfoDescription payload, independent of the array name.
    pub metadata: crate::metadata::MetaInfo,
    /// Shared processing descriptions, matching the source shared handles.
    pub data_processing: Vec<std::sync::Arc<crate::metadata::DataProcessing>>,
}

impl<T> DataArray<T> {
    pub fn new(name: impl Into<String>, data: Vec<T>) -> Self {
        Self {
            name: name.into(),
            data,
            metadata: Default::default(),
            data_processing: Vec::new(),
        }
    }
}

/// Inclusive minimum and maximum values for a nonempty dimension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NumericRange {
    pub min: f64,
    pub max: f64,
}

/// On-demand bounds for a spectrum. Empty dimensions have no range.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpectrumRanges {
    pub mz: Option<NumericRange>,
    pub intensity: Option<NumericRange>,
}

/// On-demand bounds for a chromatogram.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ChromatogramRanges {
    pub rt: Option<NumericRange>,
    pub intensity: Option<NumericRange>,
}

/// RT, m/z and intensity bounds. Each query documents which experiment data it includes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ExperimentRanges {
    pub rt: Option<NumericRange>,
    pub mz: Option<NumericRange>,
    pub intensity: Option<NumericRange>,
}

/// Spectrum representation. No automatic profile/centroid inference is made.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SpectrumType {
    #[default]
    Unknown,
    Centroid,
    Profile,
}

/// Spectrum with acquisition metadata and optional aligned annotations.
#[derive(Clone, Debug, PartialEq)]
pub struct MSSpectrum {
    pub peaks: Vec<Peak1D>,
    /// Retention time in seconds; -1 is the OpenMS default for an unset value.
    pub rt: f64,
    /// Positive MS level (1 for survey scans, 2 for MS/MS).
    pub ms_level: u32,
    pub native_id: String,
    pub name: String,
    pub spectrum_type: SpectrumType,
    pub instrument_settings: crate::metadata::InstrumentSettings,
    pub acquisition_info: crate::metadata::AcquisitionInfo,
    pub source_file: crate::metadata::SourceFile,
    /// Shared source processing handles. Conversion creates new outputs without history.
    pub data_processing: Vec<std::sync::Arc<crate::metadata::DataProcessing>>,
    pub products: Vec<crate::metadata::Product>,
    pub precursors: Vec<Precursor>,
    pub peptide_identifications: Vec<crate::identification::PeptideIdentification>,
    pub metadata: BTreeMap<String, String>,
    pub float_data_arrays: Vec<DataArray<f32>>,
    pub integer_data_arrays: Vec<DataArray<i32>>,
    pub string_data_arrays: Vec<DataArray<String>>,
}

impl Default for MSSpectrum {
    fn default() -> Self {
        Self {
            peaks: Vec::new(),
            rt: -1.0,
            ms_level: 1,
            native_id: String::new(),
            name: String::new(),
            spectrum_type: SpectrumType::Unknown,
            instrument_settings: Default::default(),
            acquisition_info: Default::default(),
            source_file: Default::default(),
            data_processing: Vec::new(),
            products: Vec::new(),
            precursors: Vec::new(),
            peptide_identifications: Vec::new(),
            metadata: BTreeMap::new(),
            float_data_arrays: Vec::new(),
            integer_data_arrays: Vec::new(),
            string_data_arrays: Vec::new(),
        }
    }
}

/// Chromatogram with retention times in seconds and aligned annotations.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MSChromatogram {
    pub instrument_settings: crate::metadata::InstrumentSettings,
    pub acquisition_info: crate::metadata::AcquisitionInfo,
    pub source_file: crate::metadata::SourceFile,
    /// Shared source processing handles. Conversion creates new outputs without history.
    pub data_processing: Vec<std::sync::Arc<crate::metadata::DataProcessing>>,
    pub chromatogram_type: crate::metadata::ChromatogramType,
    pub peaks: Vec<ChromatogramPeak>,
    pub native_id: String,
    pub name: String,
    pub precursor: Precursor,
    /// Product isolation information; XIC extraction sets its target m/z.
    pub product: crate::metadata::Product,
    pub metadata: BTreeMap<String, String>,
    pub float_data_arrays: Vec<DataArray<f32>>,
    pub integer_data_arrays: Vec<DataArray<i32>>,
    pub string_data_arrays: Vec<DataArray<String>>,
}

fn finite(value: f64, name: &str) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidValue(format!("{name} must be finite")))
    }
}

fn array_sizes<T>(arrays: &[DataArray<T>], size: usize) -> Result<()> {
    for array in arrays {
        if !array.data.is_empty() && array.data.len() != size {
            return Err(Error::InvalidValue(format!(
                "data array '{}' has {} entries for {size} peaks",
                array.name,
                array.data.len()
            )));
        }
    }
    Ok(())
}

fn select_arrays<T: Clone>(arrays: &mut [DataArray<T>], indices: &[usize]) {
    for array in arrays {
        if !array.data.is_empty() {
            array.data = indices.iter().map(|&i| array.data[i].clone()).collect();
        }
    }
}

fn validate_indices(indices: &[usize], size: usize) -> Result<()> {
    let mut seen = vec![false; size];
    for &i in indices {
        if i >= size {
            return Err(Error::InvalidValue(format!(
                "peak index {i} outside length {size}"
            )));
        }
        if seen[i] {
            return Err(Error::InvalidValue(format!("duplicate peak index {i}")));
        }
        seen[i] = true;
    }
    Ok(())
}

fn range(values: impl Iterator<Item = f64>) -> Option<NumericRange> {
    values.fold(None, |range, value| {
        Some(match range {
            None => NumericRange {
                min: value,
                max: value,
            },
            Some(NumericRange { min, max }) => NumericRange {
                min: min.min(value),
                max: max.max(value),
            },
        })
    })
}

fn check_sorted<T>(values: &[T], coordinate: impl Fn(&T) -> f64) -> Result<()> {
    for value in values {
        finite(coordinate(value), "coordinate")?;
    }
    if values
        .windows(2)
        .any(|pair| coordinate(&pair[0]) > coordinate(&pair[1]))
    {
        return Err(Error::UnsortedData);
    }
    Ok(())
}

// Searches use the first equal value. Midpoint ties select the lower position
// for peaks, while experiment RT searches select the higher position upstream.
pub(crate) fn nearest<T>(
    values: &[T],
    query: f64,
    coordinate: impl Fn(&T) -> f64,
) -> Option<usize> {
    if values.is_empty() {
        return None;
    }
    let above = values.partition_point(|value| coordinate(value) < query);
    if above == 0 {
        return Some(0);
    }
    if above == values.len() {
        return Some(above - 1);
    }
    if (coordinate(&values[above]) - query).abs() < (coordinate(&values[above - 1]) - query).abs() {
        Some(above)
    } else {
        Some(above - 1)
    }
}

// Both containers have the same rules for copying and permuting peak arrays.
macro_rules! peak_container {
    ($container:ty, $peak:ty, $position:ident) => {
        impl $container {
            pub fn new() -> Self {
                Self::default()
            }

            /// Construct a container without sorting. Validation is explicit.
            pub fn from_peaks(peaks: Vec<$peak>) -> Self {
                Self {
                    peaks,
                    ..Self::default()
                }
            }

            pub fn len(&self) -> usize {
                self.peaks.len()
            }
            pub fn is_empty(&self) -> bool {
                self.peaks.is_empty()
            }

            /// Check that every nonempty annotation array has one entry per peak.
            pub fn validate_data_arrays(&self) -> Result<()> {
                array_sizes(&self.float_data_arrays, self.len())?;
                array_sizes(&self.integer_data_arrays, self.len())?;
                array_sizes(&self.string_data_arrays, self.len())
            }

            /// Charge complete array descriptions before copying or dropping
            /// them; values are accounted separately by the consuming algorithm.
            pub(crate) fn array_descriptions_with_budget(
                &self,
                work: &mut usize,
                bytes: &mut usize,
            ) -> Result<()> {
                let count = self
                    .float_data_arrays
                    .len()
                    .checked_add(self.integer_data_arrays.len())
                    .and_then(|n| n.checked_add(self.string_data_arrays.len()))
                    .ok_or_else(|| Error::InvalidValue("array count overflow".into()))?;
                *work = work.checked_sub(count).ok_or_else(|| {
                    Error::InvalidValue("data array description work limit exceeded".into())
                })?;
                for array in &self.float_data_arrays {
                    array.description_with_budget(work, bytes)?;
                }
                for array in &self.integer_data_arrays {
                    array.description_with_budget(work, bytes)?;
                }
                for array in &self.string_data_arrays {
                    array.description_with_budget(work, bytes)?;
                }
                Ok(())
            }

            fn validate_array_descriptions(&self) -> Result<()> {
                let (mut work, mut bytes) = (50_000_000, 256 * 1024 * 1024);
                self.array_descriptions_with_budget(&mut work, &mut bytes)?;
                for array in &self.float_data_arrays {
                    array.validate_description()?;
                }
                for array in &self.integer_data_arrays {
                    array.validate_description()?;
                }
                for array in &self.string_data_arrays {
                    array.validate_description()?;
                }
                Ok(())
            }

            /// True when finite coordinates are in nondecreasing order.
            pub fn is_sorted(&self) -> bool {
                check_sorted(&self.peaks, |peak| peak.$position).is_ok()
            }

            /// Stable sort by coordinate, moving all annotation arrays together.
            /// Invalid input is rejected before any mutation.
            pub fn sort_by_position(&mut self) -> Result<()> {
                self.validate()?;
                let mut indices: Vec<usize> = (0..self.len()).collect();
                indices.sort_by(|&a, &b| {
                    self.peaks[a]
                        .$position
                        .partial_cmp(&self.peaks[b].$position)
                        .unwrap()
                });
                self.select(&indices)
            }

            /// Stable sort by intensity; `reverse` selects descending order.
            pub fn sort_by_intensity(&mut self, reverse: bool) -> Result<()> {
                self.validate()?;
                let mut indices: Vec<usize> = (0..self.len()).collect();
                indices.sort_by(|&a, &b| {
                    let order = self.peaks[a]
                        .intensity
                        .partial_cmp(&self.peaks[b].intensity)
                        .unwrap();
                    if reverse { order.reverse() } else { order }
                });
                self.select(&indices)
            }

            /// Keep/reorder unique indices and aligned arrays, preserving metadata.
            /// Bad indices or array lengths leave the container unchanged.
            pub fn select(&mut self, indices: &[usize]) -> Result<()> {
                validate_indices(indices, self.len())?;
                self.validate_data_arrays()?;
                self.peaks = indices.iter().map(|&i| self.peaks[i]).collect();
                select_arrays(&mut self.float_data_arrays, indices);
                select_arrays(&mut self.integer_data_arrays, indices);
                select_arrays(&mut self.string_data_arrays, indices);
                Ok(())
            }

            /// Retain peaks satisfying a predicate and their aligned annotations.
            pub fn retain_peaks(&mut self, mut keep: impl FnMut(&$peak) -> bool) -> Result<()> {
                self.validate_data_arrays()?;
                let indices: Vec<_> = self
                    .peaks
                    .iter()
                    .enumerate()
                    .filter_map(|(i, peak)| keep(peak).then_some(i))
                    .collect();
                self.select(&indices)
            }

            /// First peak with maximum intensity; `None` for an empty container.
            /// Like upstream, equal-intensity ties retain the first peak.
            pub fn base_peak(&self) -> Option<&$peak> {
                self.peaks.iter().reduce(|best, peak| {
                    if peak.intensity > best.intensity {
                        peak
                    } else {
                        best
                    }
                })
            }

            /// Sum intensities in storage order using OpenMS `f32` accumulation.
            /// This does not integrate over coordinate spacing.
            pub fn calculate_tic(&self) -> f32 {
                self.peaks
                    .iter()
                    .fold(0.0_f32, |sum, peak| sum + peak.intensity)
            }

            /// Clear peaks and annotation arrays. Optionally reset all metadata.
            pub fn clear(&mut self, clear_metadata: bool) {
                if clear_metadata {
                    *self = Self::default();
                } else {
                    self.peaks.clear();
                    self.float_data_arrays.clear();
                    self.integer_data_arrays.clear();
                    self.string_data_arrays.clear();
                }
            }
        }

        impl From<Vec<$peak>> for $container {
            fn from(peaks: Vec<$peak>) -> Self {
                Self::from_peaks(peaks)
            }
        }
    };
}

peak_container!(MSSpectrum, Peak1D, mz);
peak_container!(MSChromatogram, ChromatogramPeak, rt);

impl MSSpectrum {
    /// Validate finite values, positive MS level and parallel array lengths.
    /// Signed finite intensities and coordinates are permitted, as in OpenMS.
    pub fn validate(&self) -> Result<()> {
        finite(self.rt, "spectrum retention time")?;
        if self.ms_level == 0 {
            return Err(Error::InvalidValue(
                "spectrum MS level must be positive".into(),
            ));
        }
        for peak in &self.peaks {
            finite(peak.mz, "peak m/z")?;
            finite(f64::from(peak.intensity), "peak intensity")?;
        }
        for precursor in &self.precursors {
            precursor.validate()?;
        }
        for identification in &self.peptide_identifications {
            identification.validate()?;
        }
        self.validate_data_arrays()?;
        self.validate_array_descriptions()?;
        self.validate_acquisition_settings()
    }

    /// Inclusive bounds, recomputed from current peaks.
    pub fn ranges(&self) -> Result<SpectrumRanges> {
        self.validate()?;
        Ok(SpectrumRanges {
            mz: range(self.peaks.iter().map(|peak| peak.mz)),
            intensity: range(self.peaks.iter().map(|peak| f64::from(peak.intensity))),
        })
    }

    /// Index of the first peak at or above `mz`, or `len()` past the end.
    pub fn mz_begin(&self, mz: f64) -> Result<usize> {
        finite(mz, "query m/z")?;
        check_sorted(&self.peaks, |peak| peak.mz)?;
        Ok(self.peaks.partition_point(|peak| peak.mz < mz))
    }

    /// Index of the first peak strictly above `mz`, or `len()` past the end.
    pub fn mz_end(&self, mz: f64) -> Result<usize> {
        finite(mz, "query m/z")?;
        check_sorted(&self.peaks, |peak| peak.mz)?;
        Ok(self.peaks.partition_point(|peak| peak.mz <= mz))
    }

    /// Nearest index; midpoint ties choose lower m/z, exact duplicates the first.
    pub fn find_nearest(&self, mz: f64) -> Result<Option<usize>> {
        finite(mz, "query m/z")?;
        check_sorted(&self.peaks, |peak| peak.mz)?;
        Ok(nearest(&self.peaks, mz, |peak| peak.mz))
    }

    /// Nearest peak within an inclusive symmetric tolerance in Th.
    pub fn find_nearest_with_tolerance(&self, mz: f64, tolerance: f64) -> Result<Option<usize>> {
        self.find_nearest_in_window(mz, tolerance, tolerance)
    }

    /// Nearest peak within inclusive left/right tolerances in Th.
    /// If the nearest peak lies outside its side's window, try the other side.
    pub fn find_nearest_in_window(&self, mz: f64, left: f64, right: f64) -> Result<Option<usize>> {
        validate_window(mz, left, right)?;
        let Some(index) = self.find_nearest(mz)? else {
            return Ok(None);
        };
        let found = self.peaks[index].mz;
        if found >= mz - left && found <= mz + right {
            return Ok(Some(index));
        }
        let other = if found < mz {
            index.checked_add(1)
        } else {
            index.checked_sub(1)
        };
        Ok(other.filter(|&i| {
            self.peaks
                .get(i)
                .is_some_and(|peak| peak.mz >= mz - left && peak.mz <= mz + right)
        }))
    }

    /// Index of the first most-intense peak within an inclusive m/z window.
    pub fn find_highest_in_window(&self, mz: f64, left: f64, right: f64) -> Result<Option<usize>> {
        validate_window(mz, left, right)?;
        self.validate()?;
        check_sorted(&self.peaks, |peak| peak.mz)?;
        let begin = self.peaks.partition_point(|peak| peak.mz < mz - left);
        let end = self.peaks.partition_point(|peak| peak.mz <= mz + right);
        Ok((begin..end).reduce(|best, i| {
            if self.peaks[i].intensity > self.peaks[best].intensity {
                i
            } else {
                best
            }
        }))
    }
}

fn validate_window(center: f64, left: f64, right: f64) -> Result<()> {
    finite(center, "window center")?;
    finite(left, "left tolerance")?;
    finite(right, "right tolerance")?;
    if left < 0.0 || right < 0.0 {
        return Err(Error::InvalidValue("tolerances must be nonnegative".into()));
    }
    Ok(())
}

impl MSChromatogram {
    pub fn validate(&self) -> Result<()> {
        for peak in &self.peaks {
            finite(peak.rt, "chromatogram retention time")?;
            finite(f64::from(peak.intensity), "chromatogram intensity")?;
        }
        self.precursor.validate()?;
        self.product.validate()?;
        self.validate_data_arrays()?;
        self.validate_array_descriptions()?;
        self.validate_acquisition_settings()
    }

    pub fn ranges(&self) -> Result<ChromatogramRanges> {
        self.validate()?;
        Ok(ChromatogramRanges {
            rt: range(self.peaks.iter().map(|peak| peak.rt)),
            intensity: range(self.peaks.iter().map(|peak| f64::from(peak.intensity))),
        })
    }

    /// Index of the first measurement at or above `rt` (seconds).
    pub fn rt_begin(&self, rt: f64) -> Result<usize> {
        finite(rt, "query retention time")?;
        check_sorted(&self.peaks, |peak| peak.rt)?;
        Ok(self.peaks.partition_point(|peak| peak.rt < rt))
    }

    /// Index of the first measurement strictly above `rt` (seconds).
    pub fn rt_end(&self, rt: f64) -> Result<usize> {
        finite(rt, "query retention time")?;
        check_sorted(&self.peaks, |peak| peak.rt)?;
        Ok(self.peaks.partition_point(|peak| peak.rt <= rt))
    }

    /// Nearest measurement; midpoint ties choose the lower retention time.
    pub fn find_nearest(&self, rt: f64) -> Result<Option<usize>> {
        finite(rt, "query retention time")?;
        check_sorted(&self.peaks, |peak| peak.rt)?;
        Ok(nearest(&self.peaks, rt, |peak| peak.rt))
    }
}

/// Collection of spectra and chromatograms. No C++ backing library is required.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MSExperiment {
    pub spectra: Vec<MSSpectrum>,
    pub chromatograms: Vec<MSChromatogram>,
    pub metadata: BTreeMap<String, String>,
}

impl MSExperiment {
    pub fn new() -> Self {
        Self::default()
    }
    /// Number of spectra (chromatograms are separate).
    pub fn len(&self) -> usize {
        self.spectra.len()
    }
    pub fn is_empty(&self) -> bool {
        self.spectra.is_empty()
    }

    /// Find the preceding parent one MS level below, following the first
    /// precursor's native-ID reference when a matching earlier scan exists.
    /// A missing reference falls back to the most recent preceding parent.
    pub fn precursor_spectrum_index(&self, index: usize) -> Result<Option<usize>> {
        let spectrum = self.spectra.get(index).ok_or_else(|| {
            Error::InvalidValue("precursor lookup spectrum index is out of bounds".into())
        })?;
        let Some(parent_level) = spectrum.ms_level.checked_sub(1).filter(|&n| n > 0) else {
            return Ok(None);
        };
        let candidates = || {
            self.spectra[..index]
                .iter()
                .enumerate()
                .rev()
                .filter(|(_, s)| s.ms_level == parent_level)
        };
        if let Some(reference) = spectrum
            .precursors
            .first()
            .and_then(|p| p.spectrum_reference.as_deref())
        {
            if let Some((i, _)) = candidates().find(|(_, s)| s.native_id == reference) {
                return Ok(Some(i));
            }
        }
        Ok(candidates().next().map(|(i, _)| i))
    }

    pub fn validate(&self) -> Result<()> {
        for spectrum in &self.spectra {
            spectrum.validate()?;
        }
        for chromatogram in &self.chromatograms {
            chromatogram.validate()?;
        }
        Ok(())
    }

    /// Spectra sorted by RT, optionally requiring sorted m/z in every spectrum.
    pub fn is_sorted(&self, check_mz: bool) -> bool {
        check_sorted(&self.spectra, |spectrum| spectrum.rt).is_ok()
            && (!check_mz || self.spectra.iter().all(MSSpectrum::is_sorted))
    }

    /// Stable RT sort, optionally sorting m/z inside each spectrum.
    /// Validation completes before any changes are made.
    pub fn sort_spectra(&mut self, sort_mz: bool) -> Result<()> {
        for spectrum in &self.spectra {
            spectrum.validate()?;
        }
        if sort_mz {
            for spectrum in &mut self.spectra {
                spectrum.sort_by_position()?;
            }
        }
        self.spectra
            .sort_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
        Ok(())
    }

    pub fn rt_begin(&self, rt: f64) -> Result<usize> {
        finite(rt, "query retention time")?;
        check_sorted(&self.spectra, |spectrum| spectrum.rt)?;
        Ok(self.spectra.partition_point(|spectrum| spectrum.rt < rt))
    }

    pub fn rt_end(&self, rt: f64) -> Result<usize> {
        finite(rt, "query retention time")?;
        check_sorted(&self.spectra, |spectrum| spectrum.rt)?;
        Ok(self.spectra.partition_point(|spectrum| spectrum.rt <= rt))
    }

    /// Nearest spectrum index at an MS level (0 means all levels).
    /// Midpoint ties choose higher RT, matching OpenMS experiment search.
    pub fn closest_spectrum_in_rt(&self, rt: f64, ms_level: u32) -> Result<Option<usize>> {
        let begin = self.rt_begin(rt)?;
        let accepts = |spectrum: &MSSpectrum| ms_level == 0 || spectrum.ms_level == ms_level;
        let above = (begin..self.len()).find(|&i| accepts(&self.spectra[i]));
        let below = (0..begin).rev().find(|&i| accepts(&self.spectra[i]));
        Ok(match (below, above) {
            (None, right) => right,
            (left, None) => left,
            (Some(left), Some(right)) => Some(
                if rt - self.spectra[left].rt < self.spectra[right].rt - rt {
                    left
                } else {
                    right
                },
            ),
        })
    }

    /// Borrow spectra within an inclusive RT range and MS level (0 means all).
    /// Requires spectra sorted by retention time; m/z order is immaterial.
    pub fn spectra_in_rt_range(
        &self,
        min_rt: f64,
        max_rt: f64,
        ms_level: u32,
    ) -> Result<Vec<&MSSpectrum>> {
        finite(min_rt, "minimum retention time")?;
        finite(max_rt, "maximum retention time")?;
        if min_rt > max_rt {
            return Err(Error::InvalidValue("minimum RT exceeds maximum RT".into()));
        }
        let begin = self.rt_begin(min_rt)?;
        let end = self
            .spectra
            .partition_point(|spectrum| spectrum.rt <= max_rt);
        Ok(self.spectra[begin..end]
            .iter()
            .filter(|spectrum| ms_level == 0 || spectrum.ms_level == ms_level)
            .collect())
    }

    /// Distinct MS levels in ascending order, computed from current spectra.
    pub fn ms_levels(&self) -> Vec<u32> {
        self.spectra
            .iter()
            .map(|spectrum| spectrum.ms_level)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    /// Recalculate TIC from spectra at `ms_level` (1 for MS1, 0 for all).
    /// Preserves scan order and duplicate RTs; no resampling is performed.
    pub fn calculate_tic(&self, ms_level: u32) -> MSChromatogram {
        MSChromatogram::from_peaks(
            self.spectra
                .iter()
                .filter(|spectrum| ms_level == 0 || spectrum.ms_level == ms_level)
                .map(|spectrum| ChromatogramPeak::new(spectrum.rt, spectrum.calculate_tic()))
                .collect(),
        )
    }

    /// Current spectrum/peak bounds at an MS level; 0 includes all levels.
    /// Spectrum RTs are included even for empty spectra.
    pub fn ranges(&self, ms_level: u32) -> Result<ExperimentRanges> {
        let spectra: Vec<_> = self
            .spectra
            .iter()
            .filter(|spectrum| ms_level == 0 || spectrum.ms_level == ms_level)
            .collect();
        for spectrum in &spectra {
            spectrum.validate()?;
        }
        Ok(ExperimentRanges {
            rt: range(spectra.iter().map(|spectrum| spectrum.rt)),
            mz: range(
                spectra
                    .iter()
                    .flat_map(|spectrum| spectrum.peaks.iter().map(|peak| peak.mz)),
            ),
            intensity: range(
                spectra.iter().flat_map(|spectrum| {
                    spectrum.peaks.iter().map(|peak| f64::from(peak.intensity))
                }),
            ),
        })
    }

    pub fn clear(&mut self, clear_metadata: bool) {
        self.spectra.clear();
        self.chromatograms.clear();
        if clear_metadata {
            self.metadata.clear();
        }
    }
}
