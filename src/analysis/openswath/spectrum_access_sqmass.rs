// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectrum access backed by an sqMass (SQLite) spectrum store.
//!
//! Ports `ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h` with
//! `SpectrumAccessSqMass.cpp`. Requires the `sqmass` feature. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).

use super::data_structures::{ChromatogramPtr, Spectrum, SpectrumMeta, SpectrumPtr};
use super::spectrum_access::{SpectrumAccess, SpectrumAccessPtr, check_delta_rt, check_id};
use super::spectrum_access_openms::SpectrumAccessOpenMSInMemory;
use crate::format::mzml_sqlite_handler::MzMLSqliteHandler;
use crate::kernel::MSSpectrum;
use crate::{Error, Result};
use std::sync::Arc;

/// Spectrum access over an sqMass file through an [`MzMLSqliteHandler`],
/// optionally restricted to a subset of its spectra.
///
/// # Subsetting
///
/// Three constructors are provided:
///
/// - [`Self::new`] exposes every spectrum of the handler;
/// - [`Self::with_indices`] exposes only the spectra with the given SQL ids;
/// - [`Self::narrowed`] narrows an existing view further; its positions are
///   interpreted against the parent's already filtered view.
///
/// Every `*_by_id` method honours the subset: id `i` is the `i`-th element of
/// the visible view, not the `i`-th row of the SQL table.
///
/// Without a subset, an id is used **as the SQL spectrum ID**, as in the
/// source; this presumes the dense IDs `0..n` that the handler writes. A file
/// with gaps in its IDs needs an explicit subset.
///
/// # Performance and concurrency
///
/// Each [`SpectrumAccess::spectrum_by_id`] or
/// [`SpectrumAccess::spectrum_meta_by_id`] call issues one handler read and
/// pays the per-query SQLite overhead; callers that need many spectra should
/// prefer [`Self::all_spectra`], which reads in one batch. The view is
/// read-only and the handler opens its own connection per read, so
/// [`SpectrumAccess::light_clone`] shares the handler and the clone is safe to
/// use from another thread.
///
/// # Chromatograms
///
/// Only [`SpectrumAccess::nr_chromatograms`] is supported;
/// [`SpectrumAccess::chromatogram_by_id`] and
/// [`SpectrumAccess::chromatogram_native_id`] return [`Error::Unsupported`]
/// (source `Exception::NotImplemented`).
///
/// # Differences from the source
///
/// - The spectrum and chromatogram counts are read once, at construction, so
///   that [`SpectrumAccess::nr_spectra`] cannot fail; the source queries the
///   file on every call.
/// - The handler is shared by the clones (it is stateless for reads) instead of
///   copied.
/// - [`SpectrumAccess::spectrum_meta_by_id`] reads metadata only; the source
///   decodes the peak arrays too and discards them. The result is the same.
#[derive(Clone, Debug)]
pub struct SpectrumAccessSqMass {
    handler: Arc<MzMLSqliteHandler>,
    sidx: Vec<i64>,
    nr_spectra: usize,
    nr_chromatograms: usize,
}

fn light_spectrum(spectrum: &MSSpectrum) -> Spectrum {
    Spectrum::from_arrays(
        spectrum.peaks.iter().map(|p| p.mz).collect(),
        spectrum
            .peaks
            .iter()
            .map(|p| f64::from(p.intensity))
            .collect(),
    )
}

fn light_meta(spectrum: &MSSpectrum) -> SpectrumMeta {
    SpectrumMeta {
        id: spectrum.native_id.clone(),
        rt: spectrum.rt,
        ms_level: i32::try_from(spectrum.ms_level).unwrap_or(i32::MAX),
        ..SpectrumMeta::default()
    }
}

impl SpectrumAccessSqMass {
    /// A view of every spectrum in the handler's file.
    ///
    /// # Errors
    ///
    /// The handler's errors when the spectrum or chromatogram count cannot be read.
    pub fn new(handler: Arc<MzMLSqliteHandler>) -> Result<Self> {
        Self::with_indices(handler, Vec::new())
    }

    /// A view of only the spectra with the given SQL ids, in the given order.
    /// An empty list falls back to "all spectra".
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a negative id; the handler's errors when the
    /// counts cannot be read.
    pub fn with_indices(handler: Arc<MzMLSqliteHandler>, indices: Vec<i64>) -> Result<Self> {
        if let Some(bad) = indices.iter().find(|&&i| i < 0) {
            return Err(Error::InvalidValue(format!(
                "sqMass spectrum index {bad} is negative"
            )));
        }
        let nr_spectra = if indices.is_empty() {
            handler.nr_spectra()?
        } else {
            indices.len()
        };
        let nr_chromatograms = handler.nr_chromatograms()?;
        Ok(Self {
            handler,
            sidx: indices,
            nr_spectra,
            nr_chromatograms,
        })
    }

    /// A narrower view of `parent`.
    ///
    /// When `parent` has a subset, `positions` index into it (position `k`
    /// selects the parent's `k`-th visible spectrum); otherwise they are SQL
    /// ids. An empty `positions` inherits the parent's view unchanged.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a position is outside the parent's subset
    /// (source `Exception::IllegalArgument`; the source misses negative
    /// positions, which `usize` rules out here).
    pub fn narrowed(parent: &SpectrumAccessSqMass, positions: &[usize]) -> Result<Self> {
        let sidx = if positions.is_empty() {
            parent.sidx.clone()
        } else if parent.sidx.is_empty() {
            positions
                .iter()
                .map(|&p| {
                    i64::try_from(p).map_err(|_| {
                        Error::InvalidValue(format!("sqMass spectrum index {p} is too large"))
                    })
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            positions
                .iter()
                .map(|&p| {
                    parent.sidx.get(p).copied().ok_or_else(|| {
                        Error::InvalidValue(format!(
                            "error creating SpectrumAccessSqMass with an index {p} that exceeds the number of available data {}",
                            parent.sidx.len()
                        ))
                    })
                })
                .collect::<Result<Vec<_>>>()?
        };
        let nr_spectra = if sidx.is_empty() {
            parent.nr_spectra
        } else {
            sidx.len()
        };
        Ok(Self {
            handler: Arc::clone(&parent.handler),
            sidx,
            nr_spectra,
            nr_chromatograms: parent.nr_chromatograms,
        })
    }

    /// The SQL ids of the visible subset; empty when every spectrum is visible.
    pub fn subset(&self) -> &[i64] {
        &self.sidx
    }

    fn sql_id(&self, id: usize) -> Result<i64> {
        check_id(id, self.nr_spectra, "spectrum")?;
        if self.sidx.is_empty() {
            i64::try_from(id)
                .map_err(|_| Error::InvalidValue(format!("spectrum id {id} is too large")))
        } else {
            Ok(self.sidx[id])
        }
    }

    fn read_one(&self, id: usize, meta_only: bool) -> Result<MSSpectrum> {
        let sql_id = self.sql_id(id)?;
        self.handler
            .read_spectra(&[sql_id], meta_only)?
            .into_iter()
            .next()
            .ok_or_else(|| {
                Error::MissingInformation(format!("sqMass file has no spectrum with ID {sql_id}"))
            })
    }

    /// Batch-read every visible spectrum with its metadata (native ID, RT, MS
    /// level), in one read. Source `getAllSpectra`.
    ///
    /// Without a subset the whole experiment is read; with one, the subset's
    /// spectra are read in increasing SQL-ID order (the handler's order, as in
    /// the source), which differs from the subset's own order when that is not
    /// ascending.
    ///
    /// # Errors
    ///
    /// The handler's errors; a subset with a repeated id is refused by the
    /// handler.
    pub fn all_spectra(&self) -> Result<(Vec<SpectrumPtr>, Vec<SpectrumMeta>)> {
        let spectra = if self.sidx.is_empty() {
            self.handler.read_experiment(false)?.spectra
        } else {
            self.handler.read_spectra(&self.sidx, false)?
        };
        let meta = spectra.iter().map(light_meta).collect();
        let light = spectra
            .iter()
            .map(|s| Arc::new(light_spectrum(s)))
            .collect();
        Ok((light, meta))
    }

    /// A materialised copy of the visible spectra, read in one batch through
    /// [`Self::all_spectra`], with no chromatograms.
    ///
    /// This is the source's special case in the `SpectrumAccessOpenMSInMemory`
    /// constructor for an sqMass origin, made an explicit choice.
    ///
    /// # Errors
    ///
    /// Those of [`Self::all_spectra`].
    pub fn to_in_memory(&self) -> Result<SpectrumAccessOpenMSInMemory> {
        let (spectra, meta) = self.all_spectra()?;
        SpectrumAccessOpenMSInMemory::from_parts(spectra, meta, Vec::new(), Vec::new())
    }
}

impl SpectrumAccess for SpectrumAccessSqMass {
    fn light_clone(&self) -> Result<SpectrumAccessPtr> {
        Ok(Arc::new(self.clone()))
    }

    /// The m/z and intensity arrays of the spectrum; data arrays are not read.
    fn spectrum_by_id(&self, id: usize) -> Result<SpectrumPtr> {
        Ok(Arc::new(light_spectrum(&self.read_one(id, false)?)))
    }

    /// Visible positions of the spectra the handler finds for `rt ± delta_rt`.
    ///
    /// For a positive `delta_rt` the window is inclusive on both sides; for
    /// `delta_rt == 0` the handler returns at most the earliest spectrum at or
    /// after `rt` (see [`MzMLSqliteHandler::spectra_indices_by_rt`]). With a
    /// subset, each SQL id found is mapped back to every position that holds it.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for a negative or NaN `delta_rt`, and the
    /// handler's errors (which include non-finite query values).
    fn spectra_by_rt(&self, rt: f64, delta_rt: f64) -> Result<Vec<usize>> {
        check_delta_rt(rt, delta_rt)?;
        let found = self
            .handler
            .spectra_indices_by_rt(rt, delta_rt, &self.sidx)?;
        if self.sidx.is_empty() {
            return found
                .into_iter()
                .map(|id| {
                    usize::try_from(id).map_err(|_| {
                        Error::InvalidValue(format!("negative sqMass spectrum ID {id}"))
                    })
                })
                .collect();
        }
        let mut mapped = Vec::new();
        for id in found {
            for (position, &sql_id) in self.sidx.iter().enumerate() {
                if sql_id == id {
                    mapped.push(position);
                }
            }
        }
        Ok(mapped)
    }

    fn nr_spectra(&self) -> usize {
        self.nr_spectra
    }

    /// Native ID, RT and MS level of the spectrum; `index` stays 0, as in the source.
    fn spectrum_meta_by_id(&self, id: usize) -> Result<SpectrumMeta> {
        Ok(light_meta(&self.read_one(id, true)?))
    }

    /// Always [`Error::Unsupported`]: sqMass-backed access does not implement
    /// per-chromatogram payload retrieval.
    fn chromatogram_by_id(&self, _id: usize) -> Result<ChromatogramPtr> {
        Err(Error::Unsupported(
            "SpectrumAccessSqMass does not serve chromatograms".into(),
        ))
    }

    /// Total chromatogram count of the file (chromatogram subsets are not
    /// supported, as in the source).
    fn nr_chromatograms(&self) -> usize {
        self.nr_chromatograms
    }

    /// Always [`Error::Unsupported`]: sqMass-backed access does not expose
    /// chromatogram native IDs.
    fn chromatogram_native_id(&self, _id: usize) -> Result<String> {
        Err(Error::Unsupported(
            "SpectrumAccessSqMass does not serve chromatogram native IDs".into(),
        ))
    }
}
