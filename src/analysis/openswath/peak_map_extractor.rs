// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Extraction of raw m/z–RT–ion-mobility peak clouds at targeted coordinates.
//!
//! Ports `ANALYSIS/OPENSWATH/PeakMapExtractor.h` with `PeakMapExtractor.cpp`.
//! The extractor follows the targeting contract of
//! [`ChromatogramExtractorAlgorithm`](crate::analysis::openswath::chromatogram_extractor_algorithm::ChromatogramExtractorAlgorithm):
//! it takes extraction coordinates and applies their RT, m/z and optionally
//! ion-mobility windows. Unlike chromatogram extraction it keeps every matching
//! raw point instead of integrating one intensity per spectrum. The output is
//! one [`ExtractedPeakMap`](crate::analysis::openswath::peak_map_extractor::ExtractedPeakMap) per coordinate, with parallel m/z, RT, ion-mobility
//! and intensity arrays. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_EXTRACTION_SUPPORT.md).
//!
//! Unlike the chromatogram extractor's open top-hat window, the windows here
//! are **closed**: a point exactly on an edge is kept, as in the source
//! (`PeakMapExtractor.cpp:114-117`). Each coordinate is located by binary
//! search, so the coordinates need not be sorted, but each spectrum's m/z
//! array must be.

use super::chromatogram_extractor_algorithm::ExtractionCoordinates;
use super::spectrum_access::SpectrumAccess;
use crate::concept::progress_logger::{ProgressLogger, ProgressReporter, progress_value};
use crate::{Error, Result};

/// Largest total number of points one extraction may return.
///
/// Native ceiling; the source has none. Each point costs 32 bytes of output,
/// so the ceiling is 32 GiB.
pub const MAX_PEAK_MAP_POINTS: usize = 1 << 30;

/// The progress label of [`PeakMapExtractor::extract_peak_maps`], source
/// `PeakMapExtractor.cpp:50`.
pub const PEAK_MAP_PROGRESS_LABEL: &str = "Extracting targeted peak maps";

/// One extracted targeted peak map.
///
/// Source `PeakMapExtractor::ExtractedPeakMap`. The four point arrays are
/// parallel: point `i` is (`mz[i]`, `rt[i]`, `ion_mobility[i]`,
/// `intensity[i]`).
#[derive(Clone, Debug, PartialEq)]
pub struct ExtractedPeakMap {
    /// Id of the coordinate the map was extracted for.
    pub native_id: String,
    /// Target m/z of the coordinate. Default 0.
    pub target_mz: f64,
    /// Centre of the coordinate's RT window, or NaN (the default) when the
    /// window is not positive, i.e. the whole RT range was extracted.
    pub target_rt: f64,
    /// Target ion mobility of the coordinate. Default -1.
    pub target_ion_mobility: f64,
    /// RT window start of the coordinate. Default 0.
    pub rt_start: f64,
    /// RT window end of the coordinate. Default -1.
    pub rt_end: f64,
    /// m/z of each extracted point.
    pub mz: Vec<f64>,
    /// RT of the spectrum each point came from.
    pub rt: Vec<f64>,
    /// Raw ion mobility of each point, recorded even when not filtered on.
    pub ion_mobility: Vec<f64>,
    /// Intensity of each point.
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

impl ExtractedPeakMap {
    /// The empty map of a coordinate, with the target fields set as the
    /// source sets them before extracting (`PeakMapExtractor.cpp:30-41`).
    pub fn for_coordinate(coord: &ExtractionCoordinates) -> Self {
        let mut map = Self {
            native_id: coord.id.clone(),
            target_mz: coord.mz,
            target_ion_mobility: coord.ion_mobility,
            rt_start: coord.rt_start,
            rt_end: coord.rt_end,
            ..Self::default()
        };
        if map.rt_end - map.rt_start > 0.0 {
            map.target_rt = (map.rt_start + map.rt_end) / 2.0;
        }
        map
    }

    /// Number of extracted points.
    pub fn len(&self) -> usize {
        self.mz.len()
    }

    /// Whether no point was extracted.
    pub fn is_empty(&self) -> bool {
        self.mz.is_empty()
    }
}

/// Raw peak-map extraction at given coordinates.
///
/// Source `PeakMapExtractor`, a `ProgressLogger` with no state of its own; the
/// progress logger is passed to
/// [`PeakMapExtractor::extract_peak_maps_with_progress`] instead.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PeakMapExtractor;

impl PeakMapExtractor {
    /// A new extractor. Source default constructor.
    pub fn new() -> Self {
        Self
    }

    /// Extract the raw points around each coordinate.
    ///
    /// Source `extractPeakMaps`, whose output vector becomes the return value:
    /// one map per coordinate, in coordinate order, each first filled by
    /// [`ExtractedPeakMap::for_coordinate`]. A point of a spectrum is kept for
    /// a coordinate when the spectrum RT is in `[rt_start, rt_end]` (or the RT
    /// window is not positive), its m/z is in `[mz - w/2, mz + w/2]` with `w`
    /// the full `mz_extraction_window` in Th (or ppm of `mz` when `ppm`), and,
    /// when `im_extraction_window` is positive and the coordinate's
    /// `ion_mobility` is non-negative, its ion mobility is within
    /// `ion_mobility ± im_extraction_window / 2`. All windows are closed. A
    /// negative `im_extraction_window` extracts the full ion-mobility range
    /// while still recording the raw ion-mobility values. `filter` must be
    /// `"tophat"` (the source's default argument); raw points are kept, so no
    /// weighted profile applies.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when `filter` is not `"tophat"`, when a
    ///   non-empty spectrum has no ion-mobility array, or when a spectrum's
    ///   m/z, intensity and ion-mobility arrays differ in length (all the
    ///   source's `IllegalArgument`). As in the source, the filter is checked
    ///   only when there is at least one spectrum and one coordinate.
    /// - [`Error::InvalidRange`] when more than [`MAX_PEAK_MAP_POINTS`] points
    ///   would be returned.
    /// - The errors of the spectrum access.
    ///
    /// Nothing is returned on error; the source leaves its output cleared,
    /// resized and partly filled.
    pub fn extract_peak_maps(
        &self,
        input: &dyn SpectrumAccess,
        extraction_coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
    ) -> Result<Vec<ExtractedPeakMap>> {
        self.extract_reporting(
            input,
            extraction_coordinates,
            mz_extraction_window,
            ppm,
            im_extraction_window,
            filter,
            &mut ProgressReporter::silent(),
        )
    }

    /// [`PeakMapExtractor::extract_peak_maps`], reporting to `progress` as the
    /// source's `ProgressLogger` base does: `startProgress(0, n, "Extracting
    /// targeted peak maps")`, `setProgress(i)` before spectrum `i`,
    /// `endProgress()`.
    ///
    /// # Errors
    ///
    /// As [`PeakMapExtractor::extract_peak_maps`], and the errors of
    /// `progress`.
    #[allow(clippy::too_many_arguments)]
    pub fn extract_peak_maps_with_progress(
        &self,
        input: &dyn SpectrumAccess,
        extraction_coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
        progress: &mut ProgressLogger,
    ) -> Result<Vec<ExtractedPeakMap>> {
        self.extract_reporting(
            input,
            extraction_coordinates,
            mz_extraction_window,
            ppm,
            im_extraction_window,
            filter,
            &mut ProgressReporter::new(Some(progress)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn extract_reporting(
        &self,
        input: &dyn SpectrumAccess,
        coordinates: &[ExtractionCoordinates],
        mz_extraction_window: f64,
        ppm: bool,
        im_extraction_window: f64,
        filter: &str,
        reporter: &mut ProgressReporter<'_>,
    ) -> Result<Vec<ExtractedPeakMap>> {
        let input_size = input.nr_spectra();
        let mut output: Vec<ExtractedPeakMap> = coordinates
            .iter()
            .map(ExtractedPeakMap::for_coordinate)
            .collect();
        if input_size < 1 || coordinates.is_empty() {
            return Ok(output);
        }
        if filter != "tophat" {
            return Err(Error::InvalidValue("Filter needs to be tophat".into()));
        }
        let mut total = 0usize;
        let records = progress_value(input_size)?;
        reporter.section(0, records, PEAK_MAP_PROGRESS_LABEL, |reporter| {
            for scan_idx in 0..input_size {
                reporter.set_count(scan_idx)?;
                let spectrum = input.spectrum_by_id(scan_idx)?;
                let meta = input.spectrum_meta_by_id(scan_idx)?;
                let mz_data = &spectrum.mz_array().data;
                if mz_data.is_empty() {
                    continue;
                }
                let im_data = match spectrum.drift_time_array() {
                    Some(array) => &array.data,
                    None => {
                        return Err(Error::InvalidValue(
                            "Peak map extraction requires ion mobility data on every input spectrum."
                                .into(),
                        ));
                    }
                };
                let int_data = &spectrum.intensity_array().data;
                if mz_data.len() != int_data.len() || mz_data.len() != im_data.len() {
                    return Err(Error::InvalidValue(
                        "m/z, intensity, and ion mobility arrays need to have the same length."
                            .into(),
                    ));
                }
                let current_rt = meta.rt;
                for (coord, peak_map) in coordinates.iter().zip(output.iter_mut()) {
                    if coord.rt_end - coord.rt_start > 0.0
                        && (current_rt < coord.rt_start || current_rt > coord.rt_end)
                    {
                        continue;
                    }
                    let (left, right) = if ppm {
                        (
                            coord.mz - coord.mz * mz_extraction_window / 2.0 * 1.0e-6,
                            coord.mz + coord.mz * mz_extraction_window / 2.0 * 1.0e-6,
                        )
                    } else {
                        (
                            coord.mz - mz_extraction_window / 2.0,
                            coord.mz + mz_extraction_window / 2.0,
                        )
                    };
                    let use_im_filter = im_extraction_window > 0.0 && coord.ion_mobility >= 0.0;
                    let left_im = coord.ion_mobility - im_extraction_window / 2.0;
                    let right_im = coord.ion_mobility + im_extraction_window / 2.0;

                    // :109 std::lower_bound(mz, left).
                    let mut idx = mz_data.partition_point(|&m| m < left);
                    while idx < mz_data.len() && mz_data[idx] <= right {
                        let point_mz = mz_data[idx];
                        if point_mz >= left && point_mz <= right {
                            let point_im = im_data[idx];
                            if !use_im_filter || (point_im >= left_im && point_im <= right_im) {
                                total += 1;
                                if total > MAX_PEAK_MAP_POINTS {
                                    return Err(Error::InvalidRange(format!(
                                        "peak-map extraction exceeds {MAX_PEAK_MAP_POINTS} points"
                                    )));
                                }
                                peak_map.mz.push(point_mz);
                                peak_map.rt.push(current_rt);
                                peak_map.ion_mobility.push(point_im);
                                peak_map.intensity.push(int_data[idx]);
                            }
                        }
                        idx += 1;
                    }
                }
            }
            Ok(())
        })?;
        Ok(output)
    }
}
