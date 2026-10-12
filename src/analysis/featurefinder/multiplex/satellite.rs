// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Satellite peaks and satellite data points of a filtered multiplex peak.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexSatelliteCentroided.h` and
//! `FEATUREFINDER/MultiplexSatelliteProfile.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

/// A single satellite peak of a centroided experiment (source
/// `MultiplexSatelliteCentroided`).
///
/// The satellite is part of a centroided `MSExperiment`, so the indices of its
/// spectrum and of the peak within that spectrum are sufficient to specify its
/// RT, m/z and intensity. See
/// [`MultiplexFilteredPeak`](super::filtered_peak::MultiplexFilteredPeak) and
/// [`MultiplexSatelliteProfile`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MultiplexSatelliteCentroided {
    /// Index of the spectrum in the centroided experiment (source `rt_idx_`,
    /// read by `getRTidx`).
    pub rt_idx: usize,
    /// Index of the peak within the spectrum (source `mz_idx_`, read by
    /// `getMZidx`).
    pub mz_idx: usize,
}

impl MultiplexSatelliteCentroided {
    /// The satellite at spectrum `rt_idx`, peak `mz_idx` (source constructor).
    pub const fn new(rt_idx: usize, mz_idx: usize) -> Self {
        Self { rt_idx, mz_idx }
    }
}

/// A single satellite data point of spline-interpolated profile data (source
/// `MultiplexSatelliteProfile`).
///
/// The data point is interpolated, so the triplet of RT, m/z and intensity is
/// stored explicitly. RT and intensity are `f32`, as in the source. See
/// [`MultiplexFilteredPeak`](super::filtered_peak::MultiplexFilteredPeak) and
/// [`MultiplexSatelliteCentroided`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MultiplexSatelliteProfile {
    /// Retention time of the data point (source `rt_`, read by `getRT`).
    pub rt: f32,
    /// m/z of the data point (source `mz_`, read by `getMZ`).
    pub mz: f64,
    /// Interpolated intensity (source `intensity_`, read by `getIntensity`).
    pub intensity: f32,
}

impl MultiplexSatelliteProfile {
    /// The data point at `rt`, `mz` with `intensity` (source constructor).
    pub const fn new(rt: f32, mz: f64, intensity: f32) -> Self {
        Self { rt, mz, intensity }
    }
}
