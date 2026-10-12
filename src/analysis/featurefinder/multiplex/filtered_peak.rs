// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A peak that passed all multiplex filters, with its satellites.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexFilteredPeak.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::collections::BTreeMap;

use super::satellite::{MultiplexSatelliteCentroided, MultiplexSatelliteProfile};

/// A `std::multimap<size_t, T>` keyed by pattern index: each key maps to its
/// values in insertion order, and iteration visits keys in ascending order,
/// which is the source multimap's iteration order.
pub type SatelliteMultimap<T> = BTreeMap<usize, Vec<T>>;

/// A single peak that passed all filters (source `MultiplexFilteredPeak`).
///
/// Each filter result corresponds to a successful search for a particular peak
/// pattern in the centroided data. The actual m/z shifts seen in the filter
/// result may differ from the theoretical shifts listed in the peak pattern.
///
/// A filtered peak consists of a primary peak and a set of satellite peaks. The
/// primary peak lies in the mono-isotopic mass trace of the lightest peptide of
/// the multiplet. The satellites are the peaks that form the m/z shift pattern
/// relative to the primary peak within the RT band; they are the evidence on
/// which a peak may pass the filters, and the primary peak is itself one of
/// them.
///
/// In both centroid and profile mode centroided data are filtered (profile data
/// are peak-picked first), so in both modes centroided peaks make up a filtered
/// peak ([`Self::satellite_count`]). In profile mode the spline-interpolated
/// profile data points that make up these peaks are stored as well
/// ([`Self::satellite_profile_count`]).
///
/// Satellites are stored as indices, not references: the filters work on
/// temporary "white" experiments, but all indices refer to the permanent,
/// original centroided experiment.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiplexFilteredPeak {
    /// m/z of the primary peak in Th (source `mz_`), the input of the
    /// subsequent clustering step.
    mz: f64,
    /// RT of the primary peak in seconds (source `rt_`, a `float`).
    rt: f32,
    /// Index of the primary peak within its spectrum (source `mz_idx_`), used
    /// to check the blacklist.
    mz_idx: usize,
    /// Index of the primary peak's spectrum (source `rt_idx_`).
    rt_idx: usize,
    /// Pattern index (a specific mass trace) to the satellite peaks forming
    /// the pattern (source `satellites_`). Neighbouring spectra in the RT band
    /// contribute peaks of the same mass trace, hence several per key.
    satellites: SatelliteMultimap<MultiplexSatelliteCentroided>,
    /// Pattern index to the spline-interpolated satellite data points, profile
    /// data only (source `satellites_profile_`).
    satellites_profile: SatelliteMultimap<MultiplexSatelliteProfile>,
}

impl MultiplexFilteredPeak {
    /// A primary peak at `mz`, `rt`, peak `mz_idx` of spectrum `rt_idx`, with
    /// no satellites yet (source constructor).
    pub fn new(mz: f64, rt: f32, mz_idx: usize, rt_idx: usize) -> Self {
        Self {
            mz,
            rt,
            mz_idx,
            rt_idx,
            satellites: BTreeMap::new(),
            satellites_profile: BTreeMap::new(),
        }
    }

    /// m/z of the primary peak (source `getMZ`).
    pub fn mz(&self) -> f64 {
        self.mz
    }

    /// RT of the primary peak (source `getRT`).
    pub fn rt(&self) -> f32 {
        self.rt
    }

    /// Index of the primary peak within its spectrum (source `getMZidx`).
    pub fn mz_idx(&self) -> usize {
        self.mz_idx
    }

    /// Index of the primary peak's spectrum in the experiment (source
    /// `getRTidx`).
    pub fn rt_idx(&self) -> usize {
        self.rt_idx
    }

    /// Add the satellite peak `mz_idx` of spectrum `rt_idx` to mass trace
    /// `pattern_idx` (source `addSatellite(size_t, size_t, size_t)`).
    pub fn add_satellite(&mut self, rt_idx: usize, mz_idx: usize, pattern_idx: usize) {
        self.add_satellite_from(
            MultiplexSatelliteCentroided::new(rt_idx, mz_idx),
            pattern_idx,
        );
    }

    /// Add `satellite` to mass trace `pattern_idx` (source
    /// `addSatellite(const MultiplexSatelliteCentroided&, size_t)`).
    pub fn add_satellite_from(
        &mut self,
        satellite: MultiplexSatelliteCentroided,
        pattern_idx: usize,
    ) {
        self.satellites
            .entry(pattern_idx)
            .or_default()
            .push(satellite);
    }

    /// Add a satellite data point to mass trace `pattern_idx` (source
    /// `addSatelliteProfile(float, double, float, size_t)`).
    pub fn add_satellite_profile(&mut self, rt: f32, mz: f64, intensity: f32, pattern_idx: usize) {
        self.add_satellite_profile_from(
            MultiplexSatelliteProfile::new(rt, mz, intensity),
            pattern_idx,
        );
    }

    /// Add `satellite` to mass trace `pattern_idx` (source
    /// `addSatelliteProfile(const MultiplexSatelliteProfile&, size_t)`).
    pub fn add_satellite_profile_from(
        &mut self,
        satellite: MultiplexSatelliteProfile,
        pattern_idx: usize,
    ) {
        self.satellites_profile
            .entry(pattern_idx)
            .or_default()
            .push(satellite);
    }

    /// Whether the peak `mz_idx` of spectrum `rt_idx` is already a satellite, in
    /// any mass trace (source `checkSatellite`).
    pub fn check_satellite(&self, rt_idx: usize, mz_idx: usize) -> bool {
        self.satellites
            .values()
            .flatten()
            .any(|s| s.rt_idx == rt_idx && s.mz_idx == mz_idx)
    }

    /// All satellite peaks by mass trace (source `getSatellites`).
    pub fn satellites(&self) -> &SatelliteMultimap<MultiplexSatelliteCentroided> {
        &self.satellites
    }

    /// All satellite data points by mass trace (source `getSatellitesProfile`).
    pub fn satellites_profile(&self) -> &SatelliteMultimap<MultiplexSatelliteProfile> {
        &self.satellites_profile
    }

    /// The satellite peaks as `(pattern index, satellite)` pairs, in the
    /// source multimap's iteration order.
    pub fn satellite_entries(
        &self,
    ) -> impl Iterator<Item = (usize, &MultiplexSatelliteCentroided)> {
        self.satellites
            .iter()
            .flat_map(|(&key, values)| values.iter().map(move |v| (key, v)))
    }

    /// The satellite data points as `(pattern index, data point)` pairs, in the
    /// source multimap's iteration order.
    pub fn satellite_profile_entries(
        &self,
    ) -> impl Iterator<Item = (usize, &MultiplexSatelliteProfile)> {
        self.satellites_profile
            .iter()
            .flat_map(|(&key, values)| values.iter().map(move |v| (key, v)))
    }

    /// Number of satellite peaks (source `size`).
    pub fn satellite_count(&self) -> usize {
        self.satellites.values().map(Vec::len).sum()
    }

    /// Number of satellite data points (source `sizeProfile`).
    pub fn satellite_profile_count(&self) -> usize {
        self.satellites_profile.values().map(Vec::len).sum()
    }
}
