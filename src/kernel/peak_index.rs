// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{MSExperiment, MSSpectrum, Peak1D};
use crate::{Error, Result};

/// A source-compatible peak or feature index. `usize::MAX` marks an unset
/// component. Validity checks only the peak component, as in the source;
/// accessors independently check the dimensions they actually use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PeakIndex {
    pub peak: usize,
    pub spectrum: usize,
}

impl Default for PeakIndex {
    fn default() -> Self {
        Self {
            peak: usize::MAX,
            spectrum: usize::MAX,
        }
    }
}

impl PeakIndex {
    /// Index into a spectrum of a peak map; argument order follows C++.
    pub const fn new(spectrum: usize, peak: usize) -> Self {
        Self { peak, spectrum }
    }

    /// Index into a feature or consensus-feature slice; spectrum is unset.
    pub const fn for_feature(peak: usize) -> Self {
        Self {
            peak,
            spectrum: usize::MAX,
        }
    }

    /// Source validity is an index-state check, not a bounds check against a map.
    pub const fn is_valid(&self) -> bool {
        self.peak != usize::MAX
    }

    pub fn clear(&mut self) {
        *self = Self::default();
    }

    /// Checked access to any feature-like slice, ignoring the spectrum component.
    pub fn get_feature<'a, T>(&self, features: &'a [T]) -> Result<&'a T> {
        features
            .get(self.peak)
            .ok_or_else(|| Error::InvalidValue("feature index exceeds map size".into()))
    }

    /// Checked spectrum access, independent of whether the peak component is set.
    pub fn get_spectrum<'a>(&self, experiment: &'a MSExperiment) -> Result<&'a MSSpectrum> {
        experiment
            .spectra
            .get(self.spectrum)
            .ok_or_else(|| Error::InvalidValue("spectrum index exceeds map size".into()))
    }

    /// Checked access to a peak, validating spectrum before peak like the source.
    pub fn get_peak<'a>(&self, experiment: &'a MSExperiment) -> Result<&'a Peak1D> {
        self.get_spectrum(experiment)?
            .peaks
            .get(self.peak)
            .ok_or_else(|| Error::InvalidValue("peak index exceeds spectrum size".into()))
    }
}
