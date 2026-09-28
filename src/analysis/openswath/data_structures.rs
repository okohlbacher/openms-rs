// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The light data structures used by the OpenSWATH interfaces.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/DataStructures.h` (header-only). Many of these
//! types are closely related to the ProteoWizard data structures, originally
//! written by Darren Kessner and released under the Apache 2.0 licence. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md)
//! for the member mapping and evidence.
//!
//! - [`BinaryDataArray`] holds a `Vec<f64>` with the data and a description.
//! - [`ChromatogramMeta`] is the identifying information of a chromatogram.
//! - [`Chromatogram`] holds binary data arrays; the first is the time (RT)
//!   array, the second the intensity array.
//! - [`SpectrumMeta`] is the identifying information of a spectrum (index,
//!   identifier, RT, MS level).
//! - [`Spectrum`] holds binary data arrays; the first is the m/z array, the
//!   second the intensity array, and any further arrays are extra per-peak
//!   data such as ion mobility.
//!
//! # Ownership
//!
//! The source stores every array behind a `std::shared_ptr` (`BinaryDataArrayPtr`)
//! and every spectrum and chromatogram behind another (`SpectrumPtr`,
//! `ChromatogramPtr`). The port shares whole spectra and chromatograms through
//! [`SpectrumPtr`] and [`ChromatogramPtr`] (`Arc`), but a record owns its arrays
//! by value: no source consumer relies on two records aliasing one array, and
//! owning them removes the aliasing hazard in which mutating one spectrum's array
//! silently mutates another's.
//!
//! # The two leading arrays always exist
//!
//! The source getters document that the m/z (time) and intensity arrays "may be
//! null", because `setDataArrays` accepts any vector, including an empty one, and
//! `setMZArray(nullptr)` is legal. The port holds no null array: both leading
//! arrays exist from construction, and [`Spectrum::set_data_arrays`] and
//! [`Chromatogram::set_data_arrays`] refuse a list shorter than two, which is the
//! only way the source could reach a null leading array other than by storing a
//! null pointer explicitly.

use crate::{Error, Result};
use std::sync::Arc;

/// The structure into which encoded binary data goes.
///
/// The source keeps a commented-out `dataProcessingPtr` member; it has no
/// counterpart here either.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BinaryDataArray {
    /// The binary data.
    pub data: Vec<f64>,
    /// Optional data description for non-standard arrays, such as
    /// `"Ion Mobility"`; empty for the leading m/z, time and intensity arrays.
    pub description: String,
}

impl BinaryDataArray {
    /// An array holding `data` with an empty description.
    pub fn new(data: Vec<f64>) -> Self {
        Self {
            data,
            description: String::new(),
        }
    }

    /// An array holding `data` with the given description.
    pub fn with_description(data: Vec<f64>, description: impl Into<String>) -> Self {
        Self {
            data,
            description: description.into(),
        }
    }
}

/// Identifying information for a chromatogram.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ChromatogramMeta {
    /// The zero-based, consecutive index of the chromatogram in the chromatogram list.
    pub index: usize,
    /// A unique identifier for this chromatogram.
    pub id: String,
}

/// Shared handle to a [`Chromatogram`], as the source `ChromatogramPtr`.
pub type ChromatogramPtr = Arc<Chromatogram>;

/// Shared handle to a [`ChromatogramMeta`], as the source `ChromatogramMetaPtr`.
pub type ChromatogramMetaPtr = Arc<ChromatogramMeta>;

/// Shared handle to a [`Spectrum`], as the source `SpectrumPtr`.
pub type SpectrumPtr = Arc<Spectrum>;

/// Shared handle to a [`SpectrumMeta`], as the source `SpectrumMetaPtr`.
pub type SpectrumMetaPtr = Arc<SpectrumMeta>;

fn check_leading(arrays: &[BinaryDataArray], what: &str) -> Result<()> {
    if arrays.len() < 2 {
        return Err(Error::InvalidValue(format!(
            "an OpenSWATH {what} needs at least two binary data arrays, got {}",
            arrays.len()
        )));
    }
    Ok(())
}

/// A single chromatogram: a time (RT) array, an intensity array and optional
/// further arrays.
///
/// The source keeps `defaultArrayLength` (always 2) as a private member and
/// comments out precursor, product and data-processing members; none of them is
/// observable, so none is ported.
#[derive(Clone, Debug, PartialEq)]
pub struct Chromatogram {
    arrays: Vec<BinaryDataArray>,
}

impl Default for Chromatogram {
    fn default() -> Self {
        Self::new()
    }
}

impl Chromatogram {
    /// An empty chromatogram with an empty time array and an empty intensity array.
    pub fn new() -> Self {
        Self {
            arrays: vec![BinaryDataArray::default(), BinaryDataArray::default()],
        }
    }

    /// A chromatogram from a time and an intensity array.
    ///
    /// The arrays are not checked for equal length, as in the source, which
    /// never checks it either; consumers that pair them must check.
    pub fn from_arrays(time: Vec<f64>, intensity: Vec<f64>) -> Self {
        Self {
            arrays: vec![BinaryDataArray::new(time), BinaryDataArray::new(intensity)],
        }
    }

    /// The time (RT) array, always present (see the module notes on null arrays).
    pub fn time_array(&self) -> &BinaryDataArray {
        &self.arrays[0]
    }

    /// Mutable access to the time (RT) array.
    pub fn time_array_mut(&mut self) -> &mut BinaryDataArray {
        &mut self.arrays[0]
    }

    /// Replace the time array.
    pub fn set_time_array(&mut self, data: BinaryDataArray) {
        self.arrays[0] = data;
    }

    /// The intensity array, always present.
    pub fn intensity_array(&self) -> &BinaryDataArray {
        &self.arrays[1]
    }

    /// Mutable access to the intensity array.
    pub fn intensity_array_mut(&mut self) -> &mut BinaryDataArray {
        &mut self.arrays[1]
    }

    /// Replace the intensity array.
    pub fn set_intensity_array(&mut self, data: BinaryDataArray) {
        self.arrays[1] = data;
    }

    /// Read access to all underlying data arrays, the two leading ones first.
    pub fn data_arrays(&self) -> &[BinaryDataArray] {
        &self.arrays
    }

    /// Mutable access to the underlying data arrays.
    ///
    /// The source returns the vector itself; the port returns a slice so the
    /// two leading arrays cannot be removed. Append with [`Self::push_data_array`].
    pub fn data_arrays_mut(&mut self) -> &mut [BinaryDataArray] {
        &mut self.arrays
    }

    /// Append a further data array, as the source `getDataArrays().push_back(...)`.
    pub fn push_data_array(&mut self, data: BinaryDataArray) {
        self.arrays.push(data);
    }

    /// Set all binary data arrays.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `arrays` holds fewer than two arrays; the
    /// chromatogram is then unchanged. The source accepts any vector and a later
    /// `getTimeArray()` reads out of bounds.
    pub fn set_data_arrays(&mut self, arrays: Vec<BinaryDataArray>) -> Result<()> {
        check_leading(&arrays, "chromatogram")?;
        self.arrays = arrays;
        Ok(())
    }
}

/// Identifying information for a spectrum.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumMeta {
    /// The zero-based, consecutive index of the spectrum in the spectrum list.
    pub index: usize,
    /// A unique identifier for this spectrum (its native ID).
    pub id: String,
    /// Retention time in seconds.
    ///
    /// The source default constructor leaves `RT` uninitialised; the port
    /// defaults it to 0.
    pub rt: f64,
    /// MS level. The source default constructor leaves `ms_level`
    /// uninitialised; the port defaults it to 0.
    pub ms_level: i32,
}

impl SpectrumMeta {
    /// Comparator for the retention time, as the source functor `SpectrumMeta::RTLess`.
    ///
    /// Returns `a.rt < b.rt`; like the source it is not a total order when an RT
    /// is NaN.
    pub fn rt_less(a: &SpectrumMeta, b: &SpectrumMeta) -> bool {
        a.rt < b.rt
    }
}

/// Whether a data-array description names an ion-mobility array.
///
/// Ion-mobility arrays have different names depending on their source:
///
/// - `"Ion Mobility"` / `"Ion Mobility Centroid"` (OpenMS legacy, PeakPickerIM),
/// - `"mean inverse reduced ion mobility array"` (ProteoWizard diaPASEF),
/// - `"raw inverse reduced ion mobility array"` (BrukerTimsFile via CV MS:1003008),
/// - `"inverse reduced ion mobility"` (MSConvert legacy),
/// - `"mean ion mobility array"` (IMDataConverter for millisecond IM).
///
/// A description matches when it starts with `"Ion Mobility"`, or contains
/// `"inverse reduced ion mobility"` or `"ion mobility array"`. The match is
/// case-sensitive, as in the source: `"Ion mOBILITY"` does not match.
pub fn is_ion_mobility_description(description: &str) -> bool {
    description.starts_with("Ion Mobility")
        || description.contains("inverse reduced ion mobility")
        || description.contains("ion mobility array")
}

/// A single spectrum ("the structure that captures the generation of a peak
/// list"): an m/z array, an intensity array and optional further arrays.
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    arrays: Vec<BinaryDataArray>,
}

impl Default for Spectrum {
    fn default() -> Self {
        Self::new()
    }
}

impl Spectrum {
    /// An empty spectrum with an empty m/z array and an empty intensity array.
    pub fn new() -> Self {
        Self {
            arrays: vec![BinaryDataArray::default(), BinaryDataArray::default()],
        }
    }

    /// A spectrum from an m/z and an intensity array.
    ///
    /// The arrays are not checked for equal length, as in the source.
    pub fn from_arrays(mz: Vec<f64>, intensity: Vec<f64>) -> Self {
        Self {
            arrays: vec![BinaryDataArray::new(mz), BinaryDataArray::new(intensity)],
        }
    }

    /// The m/z array, always present (see the module notes on null arrays).
    pub fn mz_array(&self) -> &BinaryDataArray {
        &self.arrays[0]
    }

    /// Mutable access to the m/z array.
    pub fn mz_array_mut(&mut self) -> &mut BinaryDataArray {
        &mut self.arrays[0]
    }

    /// Replace the m/z array.
    pub fn set_mz_array(&mut self, data: BinaryDataArray) {
        self.arrays[0] = data;
    }

    /// The intensity array, always present.
    pub fn intensity_array(&self) -> &BinaryDataArray {
        &self.arrays[1]
    }

    /// Mutable access to the intensity array.
    pub fn intensity_array_mut(&mut self) -> &mut BinaryDataArray {
        &mut self.arrays[1]
    }

    /// Replace the intensity array.
    pub fn set_intensity_array(&mut self, data: BinaryDataArray) {
        self.arrays[1] = data;
    }

    /// Append `data` as the drift-time (ion-mobility) array.
    ///
    /// As in the source, the description is overwritten with `"Ion Mobility"`
    /// and the array is appended, not substituted: a spectrum that already has
    /// an ion-mobility array keeps it, and [`Self::drift_time_array`] keeps
    /// returning the first one.
    pub fn set_drift_time_array(&mut self, mut data: BinaryDataArray) {
        data.description = "Ion Mobility".to_string();
        self.arrays.push(data);
    }

    /// The first array whose description names ion mobility, or `None`.
    ///
    /// The rule is [`is_ion_mobility_description`]. The two leading arrays are
    /// searched too, as in the source, so a leading array given an ion-mobility
    /// description is returned.
    pub fn drift_time_array(&self) -> Option<&BinaryDataArray> {
        self.arrays
            .iter()
            .find(|array| is_ion_mobility_description(&array.description))
    }

    /// Index of the first ion-mobility array, by the rule of [`Self::drift_time_array`].
    pub fn drift_time_array_index(&self) -> Option<usize> {
        self.arrays
            .iter()
            .position(|array| is_ion_mobility_description(&array.description))
    }

    /// Read access to all underlying data arrays, the two leading ones first.
    pub fn data_arrays(&self) -> &[BinaryDataArray] {
        &self.arrays
    }

    /// Mutable access to the underlying data arrays.
    ///
    /// The source returns the vector itself; the port returns a slice so the
    /// two leading arrays cannot be removed. Append with [`Self::push_data_array`].
    pub fn data_arrays_mut(&mut self) -> &mut [BinaryDataArray] {
        &mut self.arrays
    }

    /// Append a further data array, as the source `getDataArrays().push_back(...)`.
    pub fn push_data_array(&mut self, data: BinaryDataArray) {
        self.arrays.push(data);
    }

    /// Set all binary data arrays.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `arrays` holds fewer than two arrays; the
    /// spectrum is then unchanged. The source accepts any vector and a later
    /// `getMZArray()` reads out of bounds.
    pub fn set_data_arrays(&mut self, arrays: Vec<BinaryDataArray>) -> Result<()> {
        check_leading(&arrays, "spectrum")?;
        self.arrays = arrays;
        Ok(())
    }

    /// Number of peaks, taken from the m/z array.
    ///
    /// Native convenience; the source has no size accessor.
    pub fn len(&self) -> usize {
        self.arrays[0].data.len()
    }

    /// Whether the m/z array is empty. Native convenience.
    pub fn is_empty(&self) -> bool {
        self.arrays[0].data.is_empty()
    }
}
