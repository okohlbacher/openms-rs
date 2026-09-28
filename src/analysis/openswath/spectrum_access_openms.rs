// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectrum access implementations over OpenMS data structures.
//!
//! Ports, from `ANALYSIS/OPENSWATH/DATAACCESS/`:
//!
//! - `SpectrumAccessOpenMS.h` → [`SpectrumAccessOpenMS`], a view of a shared
//!   [`MSExperiment`];
//! - `SpectrumAccessOpenMSInMemory.h` → [`SpectrumAccessOpenMSInMemory`], a
//!   fully materialised copy of any other [`SpectrumAccess`];
//! - `SpectrumAccessTransforming.h` → [`SpectrumAccessTransforming`], a wrapper
//!   that transforms every spectrum it serves;
//! - `SpectrumAccessQuadMZTransforming.h` → [`QuadMZTransform`] and
//!   [`SpectrumAccessQuadMZTransforming`], the quadratic m/z correction;
//! - `SimpleOpenMSSpectraAccessFactory.h` → [`SimpleOpenMSSpectraFactory`].
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md)
//! for the member mapping and evidence. `SpectrumAccessOpenMSCached.h` is not
//! ported: it needs the cached-mzML side-car reader (`FORMAT/CachedMzML.h`,
//! `FORMAT/HANDLERS/CachedMzMLHandler.h`), which the crate does not have yet.

use super::data_structures::{
    BinaryDataArray, Chromatogram, ChromatogramPtr, Spectrum, SpectrumMeta, SpectrumPtr,
};
use super::spectrum_access::{SpectrumAccess, SpectrumAccessPtr, check_delta_rt, check_id};
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::{Error, Result};
use std::sync::Arc;

/// Convert an [`MSSpectrum`] into an OpenSWATH [`Spectrum`] carrying every
/// numeric data array.
///
/// The m/z and intensity arrays come first; then one array per float data array
/// and one per integer data array, in that order, each under its name as the
/// description. String data arrays have no numeric form and are dropped, as in
/// the source (which only walks the float and integer arrays). Intensities are
/// widened from `f32`.
///
/// This is the conversion of the source `SpectrumAccessOpenMS::getSpectrumById`.
/// [`convert_to_spectrum`](super::data_access_helper::convert_to_spectrum) is
/// the leaner `OpenSwathDataAccessHelper::convertToSpectrumPtr`, which drops
/// every data array.
pub fn spectrum_with_data_arrays(spectrum: &MSSpectrum) -> Spectrum {
    let mz = spectrum.peaks.iter().map(|peak| peak.mz).collect();
    let intensity = spectrum
        .peaks
        .iter()
        .map(|peak| f64::from(peak.intensity))
        .collect();
    let mut out = Spectrum::from_arrays(mz, intensity);
    for array in &spectrum.float_data_arrays {
        out.push_data_array(BinaryDataArray::with_description(
            array.data.iter().map(|&v| f64::from(v)).collect(),
            array.name.clone(),
        ));
    }
    for array in &spectrum.integer_data_arrays {
        out.push_data_array(BinaryDataArray::with_description(
            array.data.iter().map(|&v| f64::from(v)).collect(),
            array.name.clone(),
        ));
    }
    out
}

/// Convert an [`MSChromatogram`] into an OpenSWATH [`Chromatogram`] carrying
/// every numeric data array, by the rule of [`spectrum_with_data_arrays`].
///
/// This is the conversion of the source `SpectrumAccessOpenMS::getChromatogramById`.
pub fn chromatogram_with_data_arrays(chromatogram: &MSChromatogram) -> Chromatogram {
    let rt = chromatogram.peaks.iter().map(|peak| peak.rt).collect();
    let intensity = chromatogram
        .peaks
        .iter()
        .map(|peak| f64::from(peak.intensity))
        .collect();
    let mut out = Chromatogram::from_arrays(rt, intensity);
    for array in &chromatogram.float_data_arrays {
        out.push_data_array(BinaryDataArray::with_description(
            array.data.iter().map(|&v| f64::from(v)).collect(),
            array.name.clone(),
        ));
    }
    for array in &chromatogram.integer_data_arrays {
        out.push_data_array(BinaryDataArray::with_description(
            array.data.iter().map(|&v| f64::from(v)).collect(),
            array.name.clone(),
        ));
    }
    out
}

fn sorted_by_rt(rts: impl Iterator<Item = f64>) -> bool {
    let mut previous = f64::NEG_INFINITY;
    for rt in rts {
        // A NaN RT makes the order undefined, so it counts as unsorted.
        if rt.is_nan() || previous > rt {
            return false;
        }
        previous = rt;
    }
    true
}

/// An implementation of the OpenSWATH spectrum access interface using OpenMS:
/// a view of a shared [`MSExperiment`].
///
/// Spectra and chromatograms are converted on every request by
/// [`spectrum_with_data_arrays`] and [`chromatogram_with_data_arrays`].
///
/// Cloning (the source copy constructor) and [`SpectrumAccess::light_clone`]
/// copy only the handle: afterwards both views point at the same experiment.
#[derive(Clone, Debug)]
pub struct SpectrumAccessOpenMS {
    experiment: Arc<MSExperiment>,
    sorted: bool,
}

impl SpectrumAccessOpenMS {
    /// A view of `experiment`.
    ///
    /// Whether the spectra are sorted by RT is recorded here, once, so that
    /// [`SpectrumAccess::spectra_by_rt`] can refuse an unsorted experiment
    /// without scanning it on every call; the source never checks.
    pub fn new(experiment: Arc<MSExperiment>) -> Self {
        let sorted = sorted_by_rt(experiment.spectra.iter().map(|spectrum| spectrum.rt));
        Self { experiment, sorted }
    }

    /// The underlying experiment.
    pub fn experiment(&self) -> &Arc<MSExperiment> {
        &self.experiment
    }

    /// The full record of spectrum `id`, source `getSpectraMetaInfo`.
    ///
    /// The source returns a copy of the spectrum's `SpectrumSettings`; the
    /// port's spectrum keeps those settings as fields of [`MSSpectrum`] itself,
    /// so this borrows the whole spectrum instead.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `id` is out of range.
    pub fn spectra_meta_info(&self, id: usize) -> Result<&MSSpectrum> {
        check_id(id, self.experiment.spectra.len(), "spectrum")?;
        Ok(&self.experiment.spectra[id])
    }

    /// The full record of chromatogram `id`, source `getChromatogramMetaInfo`.
    ///
    /// As [`Self::spectra_meta_info`], this borrows the whole chromatogram in
    /// place of a `ChromatogramSettings` copy.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `id` is out of range.
    pub fn chromatogram_meta_info(&self, id: usize) -> Result<&MSChromatogram> {
        check_id(id, self.experiment.chromatograms.len(), "chromatogram")?;
        Ok(&self.experiment.chromatograms[id])
    }
}

impl SpectrumAccess for SpectrumAccessOpenMS {
    fn light_clone(&self) -> Result<SpectrumAccessPtr> {
        Ok(Arc::new(self.clone()))
    }

    fn spectrum_by_id(&self, id: usize) -> Result<SpectrumPtr> {
        check_id(id, self.experiment.spectra.len(), "spectrum")?;
        Ok(Arc::new(spectrum_with_data_arrays(
            &self.experiment.spectra[id],
        )))
    }

    /// Indices of the spectra with `rt - delta_rt <= RT <= rt + delta_rt`,
    /// except that the first spectrum at or after `rt - delta_rt` is always
    /// returned (see [`SpectrumAccess::spectra_by_rt`]). The upper bound is
    /// inclusive here, unlike in the materialised view.
    fn spectra_by_rt(&self, rt: f64, delta_rt: f64) -> Result<Vec<usize>> {
        check_delta_rt(rt, delta_rt)?;
        if !self.sorted {
            return Err(Error::UnsortedData);
        }
        let spectra = &self.experiment.spectra;
        let begin = spectra.partition_point(|spectrum| spectrum.rt < rt - delta_rt);
        if begin == spectra.len() {
            return Ok(Vec::new());
        }
        let mut result = vec![begin];
        result.extend((begin + 1..spectra.len()).take_while(|&i| spectra[i].rt <= rt + delta_rt));
        Ok(result)
    }

    fn nr_spectra(&self) -> usize {
        self.experiment.spectra.len()
    }

    /// Only `rt` and `ms_level` are filled; `index` stays 0 and `id` empty, as
    /// in the source.
    fn spectrum_meta_by_id(&self, id: usize) -> Result<SpectrumMeta> {
        check_id(id, self.experiment.spectra.len(), "spectrum")?;
        let spectrum = &self.experiment.spectra[id];
        Ok(SpectrumMeta {
            rt: spectrum.rt,
            ms_level: i32::try_from(spectrum.ms_level).unwrap_or(i32::MAX),
            ..SpectrumMeta::default()
        })
    }

    fn chromatogram_by_id(&self, id: usize) -> Result<ChromatogramPtr> {
        check_id(id, self.experiment.chromatograms.len(), "chromatogram")?;
        Ok(Arc::new(chromatogram_with_data_arrays(
            &self.experiment.chromatograms[id],
        )))
    }

    fn nr_chromatograms(&self) -> usize {
        self.experiment.chromatograms.len()
    }

    fn chromatogram_native_id(&self, id: usize) -> Result<String> {
        check_id(id, self.experiment.chromatograms.len(), "chromatogram")?;
        Ok(self.experiment.chromatograms[id].native_id.clone())
    }
}

/// A spectrum access implementation that holds the complete spectrum and
/// chromatogram payload in memory.
///
/// Constructed from any other [`SpectrumAccess`] (file-backed, on-disk cached,
/// etc.): every spectrum and chromatogram is read through the source's
/// interface during construction and cached, so later accesses do not touch the
/// disk again. The view exposes exactly the same data as the source, which is no
/// longer referenced after construction.
///
/// Cloning (the source copy constructor) and [`SpectrumAccess::light_clone`]
/// share the payloads: only the per-element handles are duplicated.
///
/// # Differences from the source
///
/// The source constructor special-cases an sqMass origin and batch-loads its
/// spectra through `SpectrumAccessSqMass::getAllSpectra`, which **skips the
/// chromatograms**: the in-memory copy of an sqMass view reports zero
/// chromatograms although the origin reports its count. The port always reads
/// through the generic interface, so it never silently drops chromatograms, and
/// a chromatogram the origin cannot serve (sqMass serves none) is an error
/// instead. Use [`Self::from_spectra`] for a spectra-only copy.
#[derive(Clone, Debug, Default)]
pub struct SpectrumAccessOpenMSInMemory {
    spectra: Vec<SpectrumPtr>,
    spectra_meta: Vec<SpectrumMeta>,
    chromatograms: Vec<ChromatogramPtr>,
    chromatogram_ids: Vec<String>,
    sorted: bool,
}

impl SpectrumAccessOpenMSInMemory {
    /// Read every spectrum and chromatogram of `origin` into memory.
    ///
    /// # Errors
    ///
    /// The first error `origin` reports; nothing is kept in that case.
    pub fn new(origin: &dyn SpectrumAccess) -> Result<Self> {
        let mut out = Self::from_spectra(origin)?;
        let n = origin.nr_chromatograms();
        out.chromatograms.reserve(n);
        out.chromatogram_ids.reserve(n);
        for i in 0..n {
            out.chromatograms.push(origin.chromatogram_by_id(i)?);
            out.chromatogram_ids.push(origin.chromatogram_native_id(i)?);
        }
        Ok(out)
    }

    /// Read every spectrum of `origin` into memory and no chromatogram.
    ///
    /// This is what the source constructor does for an sqMass origin.
    ///
    /// # Errors
    ///
    /// The first error `origin` reports.
    pub fn from_spectra(origin: &dyn SpectrumAccess) -> Result<Self> {
        let n = origin.nr_spectra();
        let mut spectra = Vec::with_capacity(n);
        let mut spectra_meta = Vec::with_capacity(n);
        for i in 0..n {
            spectra.push(origin.spectrum_by_id(i)?);
            spectra_meta.push(origin.spectrum_meta_by_id(i)?);
        }
        Self::from_parts(spectra, spectra_meta, Vec::new(), Vec::new())
    }

    /// Build the view from already materialised parts.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the spectra and their metadata, or the
    /// chromatograms and their ids, differ in length (the source's two
    /// postconditions).
    pub fn from_parts(
        spectra: Vec<SpectrumPtr>,
        spectra_meta: Vec<SpectrumMeta>,
        chromatograms: Vec<ChromatogramPtr>,
        chromatogram_ids: Vec<String>,
    ) -> Result<Self> {
        if spectra.len() != spectra_meta.len() {
            return Err(Error::InvalidValue(
                "spectra and spectrum meta data need to match".into(),
            ));
        }
        if chromatograms.len() != chromatogram_ids.len() {
            return Err(Error::InvalidValue(
                "chromatograms and chromatogram ids need to match".into(),
            ));
        }
        let sorted = sorted_by_rt(spectra_meta.iter().map(|meta| meta.rt));
        Ok(Self {
            spectra,
            spectra_meta,
            chromatograms,
            chromatogram_ids,
            sorted,
        })
    }
}

impl SpectrumAccess for SpectrumAccessOpenMSInMemory {
    fn light_clone(&self) -> Result<SpectrumAccessPtr> {
        Ok(Arc::new(self.clone()))
    }

    fn spectrum_by_id(&self, id: usize) -> Result<SpectrumPtr> {
        check_id(id, self.spectra.len(), "spectrum")?;
        Ok(Arc::clone(&self.spectra[id]))
    }

    /// Indices of the cached spectra with `rt - delta_rt <= RT < rt +
    /// delta_rt`, except that the first spectrum at or after `rt - delta_rt`
    /// is always returned (see [`SpectrumAccess::spectra_by_rt`]).
    ///
    /// The upper bound is **exclusive**, as in the source code; the source
    /// documentation claims a closed window. The source also dereferences the
    /// end iterator before checking it when the window runs past the last
    /// spectrum; the port checks first.
    fn spectra_by_rt(&self, rt: f64, delta_rt: f64) -> Result<Vec<usize>> {
        check_delta_rt(rt, delta_rt)?;
        if !self.sorted {
            return Err(Error::UnsortedData);
        }
        let meta = &self.spectra_meta;
        let begin = meta.partition_point(|m| m.rt < rt - delta_rt);
        if begin == meta.len() {
            return Ok(Vec::new());
        }
        let mut result = vec![begin];
        result.extend((begin + 1..meta.len()).take_while(|&i| meta[i].rt < rt + delta_rt));
        Ok(result)
    }

    fn nr_spectra(&self) -> usize {
        self.spectra.len()
    }

    /// The metadata captured from the origin at construction, whatever fields
    /// the origin filled.
    fn spectrum_meta_by_id(&self, id: usize) -> Result<SpectrumMeta> {
        check_id(id, self.spectra_meta.len(), "spectrum")?;
        Ok(self.spectra_meta[id].clone())
    }

    fn chromatogram_by_id(&self, id: usize) -> Result<ChromatogramPtr> {
        check_id(id, self.chromatograms.len(), "chromatogram")?;
        Ok(Arc::clone(&self.chromatograms[id]))
    }

    fn nr_chromatograms(&self) -> usize {
        self.chromatograms.len()
    }

    fn chromatogram_native_id(&self, id: usize) -> Result<String> {
        check_id(id, self.chromatogram_ids.len(), "chromatogram")?;
        Ok(self.chromatogram_ids[id].clone())
    }
}

/// A transformation applied to every spectrum a [`SpectrumAccessTransforming`]
/// serves.
///
/// The source expresses a transforming view as the abstract class
/// `SpectrumAccessTransforming`, whose subclasses override `getSpectrumById`
/// and `lightClone`. The port splits the two roles: the wrapper delegates, the
/// transform transforms.
pub trait SpectrumTransform: Clone + Send + Sync + 'static {
    /// Transform `spectrum` in place.
    fn transform(&self, spectrum: &mut Spectrum);
}

/// A transforming wrapper around a spectrum access.
///
/// Every member delegates to the wrapped view, except
/// [`SpectrumAccess::spectrum_by_id`], which applies the transform to the
/// returned spectrum, and [`SpectrumAccess::light_clone`], which wraps a light
/// clone of the wrapped view with the same transform. Chromatograms and
/// spectrum metadata pass through unchanged.
///
/// # Differences from the source
///
/// The source transforms the spectrum the wrapped view returned **in place**.
/// When that view shares its spectra (the materialised
/// [`SpectrumAccessOpenMSInMemory`] returns its cached pointer), every request
/// transforms the cached spectrum again, so the second request for a spectrum
/// returns a doubly transformed m/z array and the wrapped view is corrupted.
/// The port copies a shared spectrum before transforming it, so the wrapped
/// view is never modified and every request returns the same result.
#[derive(Clone)]
pub struct SpectrumAccessTransforming<T: SpectrumTransform> {
    inner: SpectrumAccessPtr,
    transform: T,
}

impl<T: SpectrumTransform> SpectrumAccessTransforming<T> {
    /// Wrap `inner` with `transform`.
    pub fn new(inner: SpectrumAccessPtr, transform: T) -> Self {
        Self { inner, transform }
    }

    /// The wrapped view (the source's protected `sptr_`).
    pub fn inner(&self) -> &SpectrumAccessPtr {
        &self.inner
    }

    /// The transform.
    pub fn transform(&self) -> &T {
        &self.transform
    }
}

impl<T: SpectrumTransform> SpectrumAccess for SpectrumAccessTransforming<T> {
    fn light_clone(&self) -> Result<SpectrumAccessPtr> {
        Ok(Arc::new(Self {
            inner: self.inner.light_clone()?,
            transform: self.transform.clone(),
        }))
    }

    fn spectrum_by_id(&self, id: usize) -> Result<SpectrumPtr> {
        let mut spectrum = Arc::unwrap_or_clone(self.inner.spectrum_by_id(id)?);
        self.transform.transform(&mut spectrum);
        Ok(Arc::new(spectrum))
    }

    fn spectra_by_rt(&self, rt: f64, delta_rt: f64) -> Result<Vec<usize>> {
        self.inner.spectra_by_rt(rt, delta_rt)
    }

    fn nr_spectra(&self) -> usize {
        self.inner.nr_spectra()
    }

    fn spectrum_meta_by_id(&self, id: usize) -> Result<SpectrumMeta> {
        self.inner.spectrum_meta_by_id(id)
    }

    fn chromatogram_by_id(&self, id: usize) -> Result<ChromatogramPtr> {
        self.inner.chromatogram_by_id(id)
    }

    fn nr_chromatograms(&self) -> usize {
        self.inner.nr_chromatograms()
    }

    fn chromatogram_native_id(&self, id: usize) -> Result<String> {
        self.inner.chromatogram_native_id(id)
    }
}

/// A quadratic m/z transformation, the transform of the source
/// `SpectrumAccessQuadMZTransforming`.
///
/// Each m/z value is transformed with `predict = a + b * mz + c * mz^2`. This
/// can implement an on-line mass correction for TOF instruments, for example.
///
/// - With `ppm == false` the new m/z is `predict` itself.
/// - With `ppm == true`, `predict` is the predicted ppm deviation, not the new
///   mass, and the new m/z is `mz - predict * mz / 1e6`.
///
/// Only the m/z array changes; intensities and further arrays are untouched.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadMZTransform {
    /// Regression parameter 0 (constant term).
    pub a: f64,
    /// Regression parameter 1 (linear term).
    pub b: f64,
    /// Regression parameter 2 (quadratic term).
    pub c: f64,
    /// Whether the transformation is applied in the ppm domain (if false, it is
    /// applied directly in the m/z domain).
    pub ppm: bool,
}

impl QuadMZTransform {
    /// The transformed value of one m/z, by the rule on the type.
    pub fn apply(&self, mz: f64) -> f64 {
        let predict = self.a + self.b * mz + self.c * mz * mz;
        if self.ppm {
            mz - predict * mz / 1_000_000.0
        } else {
            predict
        }
    }
}

impl SpectrumTransform for QuadMZTransform {
    fn transform(&self, spectrum: &mut Spectrum) {
        for mz in &mut spectrum.mz_array_mut().data {
            *mz = self.apply(*mz);
        }
    }
}

/// A transforming m/z wrapper around spectrum access using a quadratic equation.
///
/// See [`QuadMZTransform`] for the rule and [`SpectrumAccessTransforming`] for
/// the wrapper (including where it differs from the source).
pub type SpectrumAccessQuadMZTransforming = SpectrumAccessTransforming<QuadMZTransform>;

impl SpectrumAccessQuadMZTransforming {
    /// Source constructor `SpectrumAccessQuadMZTransforming(sptr, a, b, c, ppm)`.
    ///
    /// `a`, `b` and `c` are regression parameters 0, 1 and 2; `ppm` selects
    /// whether the transformation is applied in the ppm domain (if false, it is
    /// applied directly in the m/z domain).
    pub fn quad_mz(sptr: SpectrumAccessPtr, a: f64, b: f64, c: f64, ppm: bool) -> Self {
        Self::new(sptr, QuadMZTransform { a, b, c, ppm })
    }
}

/// A factory that returns the spectrum access implementation fitting an
/// [`MSExperiment`]. Source class `SimpleOpenMSSpectraFactory`.
#[derive(Clone, Copy, Debug, Default)]
pub struct SimpleOpenMSSpectraFactory;

impl SimpleOpenMSSpectraFactory {
    /// A spectrum access for `experiment`.
    ///
    /// The source returns a `SpectrumAccessOpenMSCached` over
    /// `experiment.getLoadedFilePath()` when any spectrum or chromatogram
    /// carries a data-processing entry with the meta value `cached_data` (a
    /// metadata-only experiment loaded from a cached mzML), and a
    /// [`SpectrumAccessOpenMS`] otherwise.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`] for a cached experiment: the cached-mzML reader
    /// behind `SpectrumAccessOpenMSCached` is not ported yet. Returning a view
    /// of the metadata-only experiment instead would serve empty spectra.
    pub fn spectrum_access_openms_ptr(experiment: Arc<MSExperiment>) -> Result<SpectrumAccessPtr> {
        if Self::is_experiment_cached(&experiment) {
            return Err(Error::Unsupported(format!(
                "experiment '{}' is a cached mzML; SpectrumAccessOpenMSCached (CachedMzML) is not ported yet",
                experiment.settings.document.loaded_file_path
            )));
        }
        Ok(Arc::new(SpectrumAccessOpenMS::new(experiment)))
    }

    /// Whether any spectrum or chromatogram of `experiment` carries a
    /// data-processing entry with the meta value `cached_data`.
    ///
    /// Private in the source; public here because callers need to know which
    /// branch [`Self::spectrum_access_openms_ptr`] takes.
    pub fn is_experiment_cached(experiment: &MSExperiment) -> bool {
        let cached = |processing: &[Arc<crate::metadata::DataProcessing>]| {
            processing
                .iter()
                .any(|dp| dp.metadata.contains_key("cached_data"))
        };
        experiment
            .spectra
            .iter()
            .any(|spectrum| cached(&spectrum.data_processing))
            || experiment
                .chromatograms
                .iter()
                .any(|chromatogram| cached(&chromatogram.data_processing))
    }
}
