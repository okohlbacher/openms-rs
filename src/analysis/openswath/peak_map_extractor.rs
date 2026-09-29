// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Targeted m/z–RT–ion-mobility peak maps (`ANALYSIS/OPENSWATH/PeakMapExtractor.h`).
//!
//! Partial port: the result type [`ExtractedPeakMap`], which
//! [`XIPMParquetConsumer`](crate::format::xipm_parquet_consumer::XIPMParquetConsumer)
//! writes. The extractor itself, `PeakMapExtractor::extractPeakMaps`, is not
//! ported yet: it takes `ChromatogramExtractorAlgorithm::ExtractionCoordinates`
//! (`ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h`, unported) and
//! reads spectra through
//! [`SpectrumAccess`](crate::analysis::openswath::spectrum_access::SpectrumAccess).

/// One extracted targeted peak map.
///
/// Ports `PeakMapExtractor::ExtractedPeakMap`: the target of one extraction
/// coordinate and every raw point that fell into its window, as four parallel
/// arrays. The defaults are the source's: target RT NaN (unset), target ion
/// mobility -1 (unset), and an RT window from 0 to -1 (empty, so unset).
#[derive(Clone, Debug, PartialEq)]
pub struct ExtractedPeakMap {
    /// The coordinate's id: a transition name, or a precursor id
    /// (`<group>_Precursor_i<n>`).
    pub native_id: String,
    /// The target m/z.
    pub target_mz: f64,
    /// The target RT; NaN when unset. The extractor sets it to the window's
    /// midpoint when the window is not empty.
    pub target_rt: f64,
    /// The target ion mobility; negative when unset.
    pub target_ion_mobility: f64,
    /// The RT window's start.
    pub rt_start: f64,
    /// The RT window's end; a window whose end is not after its start is
    /// unset.
    pub rt_end: f64,
    /// m/z of every extracted point.
    pub mz: Vec<f64>,
    /// Retention time of every extracted point.
    pub rt: Vec<f64>,
    /// Ion mobility of every extracted point.
    pub ion_mobility: Vec<f64>,
    /// Intensity of every extracted point.
    pub intensity: Vec<f64>,
}

impl Default for ExtractedPeakMap {
    fn default() -> Self {
        Self {
            native_id: String::new(),
            target_mz: 0.0,
            target_rt: f64::NAN,
            target_ion_mobility: -1.0,
            rt_start: 0.0,
            rt_end: -1.0,
            mz: Vec::new(),
            rt: Vec::new(),
            ion_mobility: Vec::new(),
            intensity: Vec::new(),
        }
    }
}
