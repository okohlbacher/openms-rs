// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Number of MS spectra per MS level as a QC metric (source
//! `QC/SpectrumCount.h`). See `docs/QC_SUPPORT.md`.

use super::qc_base::{QCBase, Requires, Status};
use crate::kernel::MSExperiment;
use std::collections::BTreeMap;

/// Number of MS spectra per MS level (SpectrumCount) as a QC metric.
///
/// Source `SpectrumCount`, a stateless metric; its constructor and destructor
/// are [`Default`] and `Drop`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SpectrumCount;

impl SpectrumCount {
    /// Source `getName`: `SpectrumCount`.
    pub const NAME: &'static str = "SpectrumCount";

    /// A new metric, source `SpectrumCount()`.
    pub fn new() -> Self {
        Self
    }

    /// The number of spectra of each MS level in `exp`, keyed by MS level in
    /// ascending order, source `compute`. Levels without spectra are absent,
    /// so an empty experiment gives an empty map.
    ///
    /// The source returns `std::map<Size, UInt>`; its caller `MzQCFile` reads
    /// levels 1 and 2 with `operator[]`, which yields 0 for a missing level.
    /// A count beyond `u32::MAX` saturates here, where the source's `UInt`
    /// would wrap.
    pub fn compute(&self, exp: &MSExperiment) -> BTreeMap<usize, u32> {
        let mut counts = BTreeMap::new();
        for spectrum in &exp.spectra {
            let count: &mut u32 = counts.entry(spectrum.ms_level as usize).or_insert(0);
            *count = count.saturating_add(1);
        }
        counts
    }
}

impl QCBase for SpectrumCount {
    fn name(&self) -> &str {
        Self::NAME
    }
    /// Raw mzML data, source `requirements`.
    fn requirements(&self) -> Status {
        Status::from_flag(Requires::RawMzml)
    }
}
