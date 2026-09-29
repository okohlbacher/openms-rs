// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! One SWATH map: a spectrum view with its isolation window.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/SwathMap.h` (with its empty `SwathMap.cpp`).
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).

use super::spectrum_access::SpectrumAccessPtr;
use std::fmt;

/// Data structure holding one SWATH map with its upper and lower isolation
/// window and whether the map is MS1 or MS2.
///
/// Ion-mobility bounds of -1 mean "no ion-mobility window".
#[derive(Clone)]
pub struct SwathMap {
    /// The spectra of the map; `None` until attached. The source default
    /// constructor leaves its `shared_ptr` null.
    pub sptr: Option<SpectrumAccessPtr>,
    /// Lower isolation-window bound (m/z).
    pub lower: f64,
    /// Upper isolation-window bound (m/z).
    pub upper: f64,
    /// Isolation-window centre (m/z).
    pub center: f64,
    /// Lower ion-mobility bound; -1 when unset.
    pub im_lower: f64,
    /// Upper ion-mobility bound; -1 when unset.
    pub im_upper: f64,
    /// Whether this is the MS1 map.
    pub ms1: bool,
}

impl Default for SwathMap {
    fn default() -> Self {
        Self {
            sptr: None,
            lower: 0.0,
            upper: 0.0,
            center: 0.0,
            im_lower: -1.0,
            im_upper: -1.0,
            ms1: false,
        }
    }
}

impl fmt::Debug for SwathMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SwathMap")
            .field("sptr", &self.sptr.as_ref().map(|_| "SpectrumAccess"))
            .field("lower", &self.lower)
            .field("upper", &self.upper)
            .field("center", &self.center)
            .field("im_lower", &self.im_lower)
            .field("im_upper", &self.im_upper)
            .field("ms1", &self.ms1)
            .finish()
    }
}

impl SwathMap {
    /// A map with an m/z window and no ion-mobility window and no spectra.
    pub fn new(mz_start: f64, mz_end: f64, mz_center: f64, is_ms1: bool) -> Self {
        Self {
            lower: mz_start,
            upper: mz_end,
            center: mz_center,
            ms1: is_ms1,
            ..Self::default()
        }
    }

    /// A map with an m/z window and an ion-mobility window and no spectra.
    pub fn with_im(
        mz_start: f64,
        mz_end: f64,
        mz_center: f64,
        im_lower: f64,
        im_upper: f64,
        is_ms1: bool,
    ) -> Self {
        Self {
            lower: mz_start,
            upper: mz_end,
            center: mz_center,
            im_lower,
            im_upper,
            ms1: is_ms1,
            sptr: None,
        }
    }

    /// Whether two maps describe the same window, with the source default
    /// tolerance of 1e-6. See [`Self::is_equal_with_tolerance`].
    pub fn is_equal(&self, other: &SwathMap) -> bool {
        self.is_equal_with_tolerance(other, 1e-6)
    }

    /// Whether two maps describe the same window: every bound and the centre
    /// differ by strictly less than `tolerance`, and the MS1 flags agree.
    ///
    /// The spectra are not compared, as in the source. A NaN bound makes the
    /// maps unequal.
    pub fn is_equal_with_tolerance(&self, other: &SwathMap, tolerance: f64) -> bool {
        let close = |a: f64, b: f64| (a - b).abs() < tolerance;
        close(self.lower, other.lower)
            && close(self.upper, other.upper)
            && close(self.center, other.center)
            && close(self.im_lower, other.im_lower)
            && close(self.im_upper, other.im_upper)
            && self.ms1 == other.ms1
    }
}
