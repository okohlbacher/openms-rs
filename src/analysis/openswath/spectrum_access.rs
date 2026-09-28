// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The interface of a mass spectrometry experiment as OpenSWATH sees it.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h` with
//! `OPENSWATHALGO/DATAACCESS/ISpectrumAccess.cpp`. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md)
//! for the member mapping and evidence.
//!
//! [`SpectrumAccess`] is the trait every OpenSWATH algorithm reads spectra and
//! chromatograms through; implementations live in
//! [`spectrum_access_openms`](crate::analysis::openswath::spectrum_access_openms)
//! (an in-memory [`MSExperiment`](crate::kernel::MSExperiment), a fully
//! materialised copy, and transforming wrappers) and, behind the `sqmass`
//! feature, `spectrum_access_sqmass` (an sqMass file).
//!
//! # Concurrency
//!
//! The source's data getters are non-`const` because file-backed
//! implementations seek a shared stream, and the class relies on
//! `lightClone()` to give each thread its own copy. The port's getters take
//! `&self`, the trait requires `Send + Sync`, and an implementation that needs
//! mutable state keeps it behind its own lock. [`SpectrumAccess::light_clone`]
//! is kept: it is still the cheap way to give a thread an independent view
//! (for a file-backed view, its own connection).
//!
//! # Identifiers
//!
//! The source takes `int` ids and checks `0 <= id < size` only with
//! `OPENMS_PRECONDITION`, which is compiled out of release builds; an
//! out-of-range id is then undefined behaviour. The port takes `usize` and
//! every implementation returns [`Error::InvalidValue`] for an id outside the
//! view.

use super::data_structures::{
    BinaryDataArray, ChromatogramPtr, Spectrum, SpectrumMeta, SpectrumPtr,
};
use crate::{Error, Result};
use std::sync::Arc;

/// A vector of spectrum handles that DIA scores can operate on; allows the
/// integration of only the target regions of several neighbouring spectra.
///
/// Source `OpenMS::SpectrumSequence`, re-exported as `OpenSwath::SpectrumSequence`.
pub type SpectrumSequence = Vec<SpectrumPtr>;

/// Shared handle to a spectrum access implementation, as the source `SpectrumAccessPtr`.
pub type SpectrumAccessPtr = Arc<dyn SpectrumAccess>;

/// The interface of a mass spectrometry experiment.
///
/// Required methods mirror the source's pure virtual members; the provided
/// methods are the source's non-virtual members `getSpectrumById(id,
/// drift_start, drift_end)` and both `getMultipleSpectra` overloads.
pub trait SpectrumAccess: Send + Sync {
    /// Light clone to produce a copy for concurrent read access.
    ///
    /// Guarantees a copy of the underlying object that provides thread-safe
    /// concurrent read access to the underlying data, with minimal copy
    /// overhead. To use it, each thread calls this to produce its own copy.
    ///
    /// # Errors
    ///
    /// Implementation-specific; a file-backed view fails when it cannot open
    /// its own connection. The source cannot report such a failure.
    fn light_clone(&self) -> Result<SpectrumAccessPtr>;

    /// The spectrum at index `id`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `id` is not below [`Self::nr_spectra`]; I/O
    /// or parse errors for a file-backed view.
    fn spectrum_by_id(&self, id: usize) -> Result<SpectrumPtr>;

    /// Indices of the spectra that lie within `rt ± delta_rt`.
    ///
    /// The window rule is implementation-specific in the source, and the port
    /// keeps each rule (see the implementations): the in-memory
    /// [`MSExperiment`](crate::kernel::MSExperiment) view includes the upper
    /// bound, the materialised and cached views exclude it, and every one of
    /// them returns the first spectrum at or after `rt - delta_rt` **even when
    /// it lies beyond `rt + delta_rt`**. [`Self::multiple_spectra`] depends on
    /// that last rule, since it asks for `delta_rt = 0`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `delta_rt` is negative or not a number (a
    /// source precondition), [`Error::UnsortedData`] when the view is not sorted
    /// by RT (undefined in the source).
    fn spectra_by_rt(&self, rt: f64, delta_rt: f64) -> Result<Vec<usize>>;

    /// Number of spectra available.
    fn nr_spectra(&self) -> usize;

    /// Meta information for the spectrum at index `id`.
    ///
    /// Which fields are filled is implementation-specific in the source and
    /// kept: see each implementation.
    ///
    /// # Errors
    ///
    /// As [`Self::spectrum_by_id`].
    fn spectrum_meta_by_id(&self, id: usize) -> Result<SpectrumMeta>;

    /// The chromatogram at index `id`.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `id` is not below [`Self::nr_chromatograms`];
    /// [`Error::Unsupported`] for a view that cannot serve chromatograms.
    fn chromatogram_by_id(&self, id: usize) -> Result<ChromatogramPtr>;

    /// Number of chromatograms available.
    fn nr_chromatograms(&self) -> usize;

    /// Native ID of the chromatogram at index `id`.
    ///
    /// # Errors
    ///
    /// As [`Self::chromatogram_by_id`].
    fn chromatogram_native_id(&self, id: usize) -> Result<String>;

    /// The spectrum at index `id`, filtered to the peaks whose drift time lies
    /// in `[drift_start, drift_end]`; the result is a copy.
    ///
    /// Source overload `getSpectrumById(int id, double drift_start, double drift_end)`.
    ///
    /// # Errors
    ///
    /// Those of [`Self::spectrum_by_id`] and of [`filter_by_drift`].
    fn spectrum_by_id_in_drift_range(
        &self,
        id: usize,
        drift_start: f64,
        drift_end: f64,
    ) -> Result<SpectrumPtr> {
        let spectrum = self.spectrum_by_id(id)?;
        Ok(Arc::new(filter_by_drift(&spectrum, drift_start, drift_end)?))
    }

    /// Fetch a spectrum sequence closest to `rt`.
    ///
    /// The closest spectrum is found from `spectra_by_rt(rt, 0.0)`: its first
    /// index, or the index before it when that one's RT is strictly closer to
    /// `rt`. The sequence starts with the closest spectrum and then alternates
    /// outward, lower neighbour first, for `nr_spectra_to_fetch / 2` (integer
    /// division) steps, skipping neighbours outside the view. An even count
    /// therefore fetches one spectrum more than asked for (`4` fetches up to
    /// five); this is the source rule.
    ///
    /// Source `getMultipleSpectra(double RT, int nr_spectra_to_fetch)`. The
    /// source takes an `int`; a negative count fetched only the closest
    /// spectrum, which is what `0` and `1` do here.
    ///
    /// An empty view, or an `rt` beyond the last spectrum, yields an empty sequence.
    ///
    /// # Errors
    ///
    /// Those of [`Self::spectra_by_rt`], [`Self::spectrum_meta_by_id`] and
    /// [`Self::spectrum_by_id`].
    fn multiple_spectra(&self, rt: f64, nr_spectra_to_fetch: usize) -> Result<SpectrumSequence> {
        let Some(order) = multiple_spectra_order(self, rt, nr_spectra_to_fetch)? else {
            return Ok(Vec::new());
        };
        order.into_iter().map(|id| self.spectrum_by_id(id)).collect()
    }

    /// Fetch a spectrum sequence closest to `rt` and filter every spectrum to
    /// the drift-time window `[drift_start, drift_end]`.
    ///
    /// Selection is exactly that of [`Self::multiple_spectra`]; each spectrum is
    /// then fetched through [`Self::spectrum_by_id_in_drift_range`].
    ///
    /// # Errors
    ///
    /// Those of [`Self::multiple_spectra`] and of [`filter_by_drift`].
    fn multiple_spectra_in_drift_range(
        &self,
        rt: f64,
        nr_spectra_to_fetch: usize,
        drift_start: f64,
        drift_end: f64,
    ) -> Result<SpectrumSequence> {
        let Some(order) = multiple_spectra_order(self, rt, nr_spectra_to_fetch)? else {
            return Ok(Vec::new());
        };
        order
            .into_iter()
            .map(|id| self.spectrum_by_id_in_drift_range(id, drift_start, drift_end))
            .collect()
    }
}

// The spectrum indices `getMultipleSpectra` visits, in its order.
fn multiple_spectra_order<A: SpectrumAccess + ?Sized>(
    access: &A,
    rt: f64,
    nr_spectra_to_fetch: usize,
) -> Result<Option<Vec<usize>>> {
    let indices = access.spectra_by_rt(rt, 0.0)?;
    let Some(&first) = indices.first() else {
        return Ok(None);
    };
    let mut closest = first;
    if first != 0 {
        let before = access.spectrum_meta_by_id(first - 1)?.rt;
        let at = access.spectrum_meta_by_id(first)?.rt;
        if (before - rt).abs() < (at - rt).abs() {
            closest -= 1;
        }
    }
    let total = access.nr_spectra();
    let mut order = vec![closest];
    for i in 1..=nr_spectra_to_fetch / 2 {
        if let Some(lower) = closest.checked_sub(i) {
            order.push(lower);
        }
        if closest + i < total {
            order.push(closest + i);
        }
    }
    Ok(Some(order))
}

/// Filter a spectrum to the peaks whose drift time lies in `[drift_start,
/// drift_end]` (both bounds inclusive); the returned spectrum is a copy.
///
/// The result has three arrays: m/z, intensity, and the ion-mobility array
/// under its original description. Any further arrays of the input are not
/// carried over, as in the source. Peak order is kept. The source notes that
/// this is inefficient because the ion-mobility array is unsorted, so every
/// peak is visited.
///
/// # Errors
///
/// - [`Error::MissingInformation`] when `input` has no ion-mobility array (see
///   [`Spectrum::drift_time_array`]). The source dereferences a null pointer;
///   its own guard for this case is commented out.
/// - [`Error::InvalidValue`] when the intensity or ion-mobility array is not as
///   long as the m/z array. The source walks all three with the m/z length and
///   reads out of bounds.
pub fn filter_by_drift(input: &Spectrum, drift_start: f64, drift_end: f64) -> Result<Spectrum> {
    let im = input.drift_time_array().ok_or_else(|| {
        Error::MissingInformation(
            "cannot filter by drift time: the spectrum has no ion-mobility array".into(),
        )
    })?;
    let mz = &input.mz_array().data;
    let intensity = &input.intensity_array().data;
    if intensity.len() != mz.len() || im.data.len() != mz.len() {
        return Err(Error::InvalidValue(format!(
            "cannot filter by drift time: m/z, intensity and ion-mobility arrays have {}, {} and {} entries",
            mz.len(),
            intensity.len(),
            im.data.len()
        )));
    }
    let mut mz_out = Vec::new();
    let mut intensity_out = Vec::new();
    let mut im_out = Vec::new();
    for ((&m, &i), &d) in mz.iter().zip(intensity).zip(&im.data) {
        if drift_start <= d && drift_end >= d {
            mz_out.push(m);
            intensity_out.push(i);
            im_out.push(d);
        }
    }
    let mut output = Spectrum::from_arrays(mz_out, intensity_out);
    output.push_data_array(BinaryDataArray::with_description(
        im_out,
        im.description.clone(),
    ));
    Ok(output)
}

/// Check a spectrum or chromatogram id against the size of a view.
pub(crate) fn check_id(id: usize, size: usize, what: &str) -> Result<()> {
    if id >= size {
        return Err(Error::InvalidValue(format!(
            "{what} id {id} is out of range for {size} {what}s"
        )));
    }
    Ok(())
}

/// Check the `deltaRT >= 0` precondition of `getSpectraByRT`.
pub(crate) fn check_delta_rt(rt: f64, delta_rt: f64) -> Result<()> {
    if delta_rt.is_nan() || delta_rt < 0.0 {
        return Err(Error::InvalidValue(format!(
            "delta RT must be a non-negative number, got {delta_rt}"
        )));
    }
    if rt.is_nan() {
        return Err(Error::InvalidValue("query RT must not be NaN".into()));
    }
    Ok(())
}
