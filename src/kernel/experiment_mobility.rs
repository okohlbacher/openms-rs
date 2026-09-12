// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ion mobility and dense RT/m/z views of a whole run.
//!
//! Source: the ion-mobility and rasterization half of `KERNEL/MSExperiment.h` /
//! `KERNEL/MSExperiment.cpp` at Core SDK `bc9cc12` — the scan-mobility searches
//! `IMBegin`/`IMEnd`, the frame predicate `isIMFrame`, the two mobility-carrying
//! bulk exports `get2DPeakDataIM` and `get2DPeakDataIMPerSpectrum`, and the
//! retention-time/m/z rasterizer `rasterizeRTMZ`.
//! `docs/EXPERIMENT_MOBILITY_SUPPORT.md` lists every member of the header, its
//! counterpart and the native differences.
//!
//! Two ion-mobility representations coexist in the source and both appear here.
//! A conventional spectrum carries one drift time for the whole scan; an *IM
//! frame* annotates every peak through a float data array. The searches and the
//! frame predicate read the **scalar** drift time; the bulk exports read the
//! **per-peak array**. Neither derivation is repeated in this module: the scalar
//! sentinel and the array lookup both come from
//! [`crate::kernel::spectrum_mobility`], and the mobility *ranges* of a run come
//! from [`crate::kernel::ranges`].
//!
//! Every operation whose cost scales with the run is checked against
//! [`MSExperiment::MAX_MOBILITY_ITEMS`](crate::kernel::MSExperiment::MAX_MOBILITY_ITEMS),
//! [`MSExperiment::MAX_RASTER_PIXELS`](crate::kernel::MSExperiment::MAX_RASTER_PIXELS)
//! or an explicit
//! [`ExperimentMobilityLimits`](crate::kernel::experiment_mobility::ExperimentMobilityLimits)
//! before anything is allocated, and every export is committed atomically, so a
//! rejected call leaves both the run and the caller's output unchanged.
//!
//! The source parallelises `rasterizeRTMZ` with OpenMP: it allocates one
//! full-size `f32` accumulation buffer per thread, fills them under
//! `#pragma omp parallel for` (`MSExperiment.cpp:409`, loop body to
//! `MSExperiment.cpp:487`) and merges them into the output afterwards
//! (`MSExperiment.cpp:489-517`). This port is serial: it reproduces the
//! source's own single-threaded branch, `if (num_threads <= 1)` at
//! `MSExperiment.cpp:355`, which writes straight into the output and returns at
//! `MSExperiment.cpp:398`.
//!
//! The two source branches agree exactly only for
//! [`RasterAggregation::Max`](crate::kernel::spectrum_mobility::RasterAggregation::Max),
//! because a maximum is associative and commutative. For
//! [`RasterAggregation::Sum`](crate::kernel::spectrum_mobility::RasterAggregation::Sum)
//! they need not: the parallel branch accumulates
//! into per-thread `f32` buffers (`MSExperiment.cpp:460`) and then adds those
//! partial sums into the pixel (`MSExperiment.cpp:499`), so the `f32`
//! summation order — and with it the rounded result — depends on how the
//! spectra were distributed over threads. This port matches the
//! single-threaded branch, not the parallel one. The performance gap is stated
//! rather than closed, because the crate introduces no threads.

use super::spectrum_mobility::RasterAggregation;
use super::{AreaBounds, AreaIter, AreaOptions, MSExperiment, MSSpectrum, check_sorted, finite};
use crate::metadata::DriftTimeUnit;
use crate::{Error, Result};
use std::mem::size_of;

/// The sentinel the source writes for a peak with no ion mobility value
/// (`MSExperiment.cpp:201`, `MSExperiment.cpp:257`). It is the same value as
/// `IMTypes::DRIFTTIME_NOT_SET`, so a genuine mobility of `-1` and a missing
/// one are indistinguishable in the source output; this port preserves that.
const MOBILITY_NOT_SET: f32 = -1.0;

/// The value the source's row cursor starts at (`float t = -1.0;`,
/// `MSExperiment.cpp:181`). A first selected retention time of exactly `-1`
/// therefore opens no row.
const ROW_RT_START: f32 = -1.0;

/// Bulk peak export as four parallel `f32` columns, one entry per peak.
///
/// Ports `MSExperiment::get2DPeakDataIM`. Extends
/// [`FlatPeakData`](crate::kernel::FlatPeakData) with the peak's ion mobility,
/// so all four vectors share a length and index together. Coordinates narrow to
/// `f32` as in the source.
///
/// Each peak's mobility comes from its **own** spectrum. The source declares
/// its row cursor `float t = -1.0;` *inside* the per-peak loop
/// (`MSExperiment.cpp:239-241`), so the cursor never advances and the
/// `it.getRT() != t` guard is true for every peak except one whose retention
/// time is exactly `-1` — for that peak the array is never fetched and `-1` is
/// written even though the spectrum carries ion mobility. This port always
/// fetches, because the source path discards data that cannot be recovered, and
/// `-1` is also the "no mobility" sentinel, so a caller cannot tell the two
/// apart. Since `-1` is OpenMS's unset retention time the case is reachable.
/// No source-compatibility option reproduces it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FlatPeakDataIm {
    /// Retention time of the spectrum each peak came from.
    pub rt: Vec<f32>,
    /// Mass-to-charge of each peak.
    pub mz: Vec<f32>,
    /// Intensity of each peak.
    pub intensity: Vec<f32>,
    /// Ion mobility of each peak, or `-1` when its spectrum carries none.
    pub ion_mobility: Vec<f32>,
}

/// Bulk peak export grouped into one row per retention time, with ion mobility.
///
/// Ports `MSExperiment::get2DPeakDataIMPerSpectrum`. All four vectors index
/// together, one entry per row.
///
/// Rows are **not** spectra. The source compares the `f64` retention time
/// against its `f32` narrowing per peak (`MSExperiment.cpp:184-186`), so two
/// spectra whose retention times narrow to the same `f32` merge into one row,
/// and a spectrum whose retention time is not exactly representable in `f32`
/// can split across rows. This is the same grouping
/// [`SpectrumPeakData`](crate::kernel::SpectrumPeakData) documents, and it has
/// one extra consequence here: the ion-mobility array is fetched only when a
/// row *starts* (`MSExperiment.cpp:188`), so the peaks of a second spectrum
/// merged into an existing row are looked up in the **first** spectrum's array.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumPeakDataIm {
    /// Retention time of each row.
    pub rt: Vec<f32>,
    /// Mass-to-charge values of each row's peaks.
    pub mz: Vec<Vec<f32>>,
    /// Intensities of each row's peaks.
    pub intensity: Vec<Vec<f32>>,
    /// Ion mobility of each row's peaks, or `-1` where none is available.
    pub ion_mobility: Vec<Vec<f32>>,
}

/// Cumulative ceilings for one mobility-carrying export call.
///
/// Native bounds with no source counterpart. They cover area validation, any
/// output already present when appending, the newly produced output and the
/// scratch used to build it, and are checked before anything is written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExperimentMobilityLimits {
    /// Maximum spectra visited.
    pub max_spectra: usize,
    /// Maximum peaks visited.
    pub max_peaks: usize,
    /// Maximum points written, counting output already present when appending.
    pub max_output_points: usize,
    /// Maximum rows written, counting output already present when appending.
    pub max_output_rows: usize,
    /// Maximum weighted visits and comparisons.
    pub max_work: usize,
    /// Conservative ceiling on logical payload; an estimate, not measured.
    pub max_bytes: usize,
}
impl Default for ExperimentMobilityLimits {
    fn default() -> Self {
        Self {
            max_spectra: 1_000_000,
            max_peaks: 10_000_000,
            max_output_points: 10_000_000,
            max_output_rows: 1_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}

/// The regular grid a run is rasterized onto, in retention time and m/z.
///
/// Source `MSExperiment::rasterizeRTMZ`'s scalar parameters
/// (`MSExperiment.h:468-477`), gathered into one value so the call site names
/// each bound. The m/z axis is the slow (row) axis and the retention-time axis
/// the fast (column) axis, matching the source's row-major output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RtMzRaster {
    /// Number of bins along the retention-time axis (image width). Must be
    /// positive.
    pub rt_bins: usize,
    /// Number of bins along the m/z axis (image height). Must be positive.
    pub mz_bins: usize,
    /// Inclusive lower end of the retention-time range, in seconds.
    pub min_rt: f64,
    /// Inclusive upper end of the retention-time range; must exceed `min_rt`.
    pub max_rt: f64,
    /// Inclusive lower end of the m/z range, in Th.
    pub min_mz: f64,
    /// Inclusive upper end of the m/z range; must exceed `min_mz`.
    pub max_mz: f64,
    /// MS level of the spectra to include; matched exactly, so `0` selects only
    /// scans whose level is zero rather than all of them.
    pub ms_level: u32,
    /// How peaks sharing a pixel combine.
    pub aggregation: RasterAggregation,
}

impl RtMzRaster {
    /// A grid over the given ranges at `ms_level`, summing peaks that share a
    /// pixel as the source default does.
    pub fn new(
        rt_bins: usize,
        mz_bins: usize,
        min_rt: f64,
        max_rt: f64,
        min_mz: f64,
        max_mz: f64,
        ms_level: u32,
    ) -> Self {
        Self {
            rt_bins,
            mz_bins,
            min_rt,
            max_rt,
            min_mz,
            max_mz,
            ms_level,
            aggregation: RasterAggregation::Sum,
        }
    }

    /// The same grid with a different aggregation mode.
    pub fn with_aggregation(mut self, aggregation: RasterAggregation) -> Self {
        self.aggregation = aggregation;
        self
    }

    /// Number of pixels the grid produces.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `rt_bins * mz_bins` overflows
    /// `usize`. The source multiplies unchecked (`MSExperiment.cpp:297`) and
    /// then fills `rt_bins * mz_bins` floats through a caller-supplied pointer,
    /// so an overflowing product there silently writes out of bounds.
    pub fn pixels(&self) -> Result<usize> {
        self.rt_bins
            .checked_mul(self.mz_bins)
            .ok_or_else(|| invalid("raster pixel count overflows"))
    }

    fn validate(&self) -> Result<usize> {
        if self.rt_bins == 0 {
            return Err(invalid("number of RT bins must be positive"));
        }
        if self.mz_bins == 0 {
            return Err(invalid("number of m/z bins must be positive"));
        }
        finite(self.min_rt, "minimum retention time")?;
        finite(self.max_rt, "maximum retention time")?;
        finite(self.min_mz, "minimum m/z")?;
        finite(self.max_mz, "maximum m/z")?;
        if self.min_rt >= self.max_rt {
            return Err(Error::InvalidRange(
                "minimum retention time must be below the maximum".into(),
            ));
        }
        if self.min_mz >= self.max_mz {
            return Err(Error::InvalidRange(
                "minimum m/z must be below the maximum".into(),
            ));
        }
        let pixels = self.pixels()?;
        if pixels > MSExperiment::MAX_RASTER_PIXELS {
            return Err(invalid("raster pixel count exceeds MAX_RASTER_PIXELS"));
        }
        Ok(pixels)
    }
}

impl MSExperiment {
    /// Largest number of spectra or peaks the scan-mobility operations in this
    /// module will visit. Matches
    /// [`MSSpectrum::MAX_MOBILITY_ITEMS`](crate::kernel::MSSpectrum::MAX_MOBILITY_ITEMS).
    pub const MAX_MOBILITY_ITEMS: usize = 100_000_000;

    /// Largest raster [`Self::rasterize_rt_mz`] will allocate, in pixels. At
    /// `f32` per pixel this caps one image at 64 MiB, as for
    /// [`MSSpectrum::MAX_RASTER_PIXELS`](crate::kernel::MSSpectrum::MAX_RASTER_PIXELS).
    pub const MAX_RASTER_PIXELS: usize = 16_777_216;

    /// Index of the first spectrum whose scan-wide drift time is not less than
    /// `im`, or `len()` past the end.
    ///
    /// Ports `IMBegin` (`MSExperiment.cpp:642-647`), a `lower_bound` under
    /// `MSSpectrum::IMLess`, which compares `getDriftTime()`
    /// (`MSSpectrum.cpp:803-806`). The source declares only the const overload
    /// even though its class test names `Iterator IMBegin(CoordinateType im)`;
    /// an index serves both.
    ///
    /// Only the scalar drift time is consulted. A spectrum whose peaks carry an
    /// ion-mobility float data array but whose scalar drift time is unset sorts
    /// at the sentinel `-1`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when the spectra are not ordered by
    /// drift time, which the source's `@note` states as a precondition without
    /// checking it, and [`Error::InvalidValue`] when `im` or a stored drift
    /// time is not finite, or when the run exceeds
    /// [`Self::MAX_MOBILITY_ITEMS`] spectra.
    pub fn im_begin(&self, im: f64) -> Result<usize> {
        finite(im, "query ion mobility")?;
        self.check_mobility_order()?;
        Ok(self
            .spectra
            .partition_point(|spectrum| spectrum.drift_time < im))
    }

    /// Index of the first spectrum whose scan-wide drift time is greater than
    /// `im`, or `len()` past the end.
    ///
    /// Ports `IMEnd` (`MSExperiment.cpp:649-654`), an `upper_bound` under the
    /// same comparator.
    ///
    /// # Errors
    ///
    /// As [`Self::im_begin`].
    pub fn im_end(&self, im: f64) -> Result<usize> {
        finite(im, "query ion mobility")?;
        self.check_mobility_order()?;
        Ok(self
            .spectra
            .partition_point(|spectrum| spectrum.drift_time <= im))
    }

    /// Whether all spectra form a single ion-mobility frame: one retention time
    /// throughout, with the drift time changing from each scan to the next.
    ///
    /// Ports `isIMFrame` (`MSExperiment.cpp:1322-1333`). An empty run is not a
    /// frame.
    ///
    /// Two source details are preserved exactly. The retention time of the
    /// *first* spectrum is the reference, compared with `!=`, so a run whose
    /// scans share a retention time only approximately is not a frame. And the
    /// drift time is compared against the **immediately preceding** scan only,
    /// not against every earlier one, so drift times `1, 2, 1` still report a
    /// frame while `1, 1, 2` does not; the source never sorts or deduplicates.
    ///
    /// A single spectrum therefore always reports `true` when its retention
    /// time is finite, even with no ion mobility information at all, because
    /// the initial "previous drift time" is `f64::MIN` and any stored value —
    /// including the unset sentinel `-1` — differs from it. That is the source's
    /// answer and it is reproduced here; it is recorded as a source defect
    /// candidate in the support document.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a retention time or drift time is
    /// not finite, or when the run exceeds [`Self::MAX_MOBILITY_ITEMS`]
    /// spectra. The source compares whatever it holds, and a NaN retention time
    /// makes its first `!=` test true, so a single NaN scan is reported as not
    /// a frame rather than rejected.
    pub fn is_im_frame(&self) -> Result<bool> {
        self.cap_spectra()?;
        let Some(first) = self.spectra.first() else {
            return Ok(false);
        };
        finite(first.rt, "spectrum retention time")?;
        let mut previous = f64::MIN;
        for spectrum in &self.spectra {
            finite(spectrum.rt, "spectrum retention time")?;
            finite(spectrum.drift_time, "spectrum drift time")?;
            if spectrum.rt != first.rt {
                return Ok(false);
            }
            if spectrum.drift_time == previous {
                return Ok(false);
            }
            previous = spectrum.drift_time;
        }
        Ok(true)
    }

    /// Every peak inside `bounds` at `ms_level`, as four parallel columns
    /// including ion mobility.
    ///
    /// Ports `get2DPeakDataIM` (`MSExperiment.cpp:226-260`). Each peak's
    /// mobility comes from its own spectrum; see [`FlatPeakDataIm`] for the
    /// source's retention-time `-1` defect this does not reproduce.
    ///
    /// # Errors
    ///
    /// As [`MSExperiment::get_2d_peak_data`], plus [`Error::InvalidValue`] when
    /// a spectrum's ion-mobility array is shorter than the peak being exported.
    pub fn get_2d_peak_data_im(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
    ) -> Result<FlatPeakDataIm> {
        self.get_2d_peak_data_im_with_limits(bounds, ms_level, ExperimentMobilityLimits::default())
    }
    /// As [`Self::get_2d_peak_data_im`], with explicit resource ceilings.
    pub fn get_2d_peak_data_im_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        limits: ExperimentMobilityLimits,
    ) -> Result<FlatPeakDataIm> {
        let mut result = FlatPeakDataIm::default();
        self.append_2d_peak_data_im_with_limits(bounds, ms_level, &mut result, limits)?;
        Ok(result)
    }
    /// Append the selected peaks to existing columns, keeping their contents.
    ///
    /// The source never clears the vectors it is handed
    /// (`MSExperiment.cpp:226-260`), so appending is the source behaviour and
    /// the fresh-output form is the convenience wrapper.
    ///
    /// # Errors
    ///
    /// As [`Self::get_2d_peak_data_im`]; the ceilings count what is already
    /// there, and misaligned column lengths are rejected.
    pub fn append_2d_peak_data_im(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut FlatPeakDataIm,
    ) -> Result<()> {
        self.append_2d_peak_data_im_with_limits(
            bounds,
            ms_level,
            output,
            ExperimentMobilityLimits::default(),
        )
    }
    /// As [`Self::append_2d_peak_data_im`], with explicit resource ceilings.
    pub fn append_2d_peak_data_im_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut FlatPeakDataIm,
        limits: ExperimentMobilityLimits,
    ) -> Result<()> {
        let mut work = Work::new(limits);
        let points = selected(self, bounds, ms_level, limits, &mut work)?;
        if points.len() == 0 {
            return Ok(());
        }
        let old = output.rt.len();
        if old != output.mz.len()
            || old != output.intensity.len()
            || old != output.ion_mobility.len()
        {
            return Err(invalid("bulk peak output vectors are not aligned"));
        }
        let count = add(old, points.len())?;
        cap(count, limits.max_output_points)?;
        work.consume(mul(count, 4)?)?;
        let mut staged = FlatPeakDataIm {
            rt: copied(&output.rt, count, &mut work)?,
            mz: copied(&output.mz, count, &mut work)?,
            intensity: copied(&output.intensity, count, &mut work)?,
            ion_mobility: copied(&output.ion_mobility, count, &mut work)?,
        };
        for point in points {
            staged.rt.push(narrow(point.spectrum.rt, "RT")?);
            staged.mz.push(narrow(point.peak.mz, "m/z")?);
            staged
                .intensity
                .push(finite_intensity(point.peak.intensity)?);
            staged
                .ion_mobility
                .push(peak_mobility(point.spectrum, point.peak_index)?);
        }
        *output = staged;
        Ok(())
    }

    /// Every peak inside `bounds` at `ms_level`, grouped into rows by
    /// `f32`-narrowed retention time and carrying ion mobility; see
    /// [`SpectrumPeakDataIm`] on grouping.
    ///
    /// # Errors
    ///
    /// As [`Self::get_2d_peak_data_im`].
    pub fn get_2d_peak_data_im_per_spectrum(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
    ) -> Result<SpectrumPeakDataIm> {
        self.get_2d_peak_data_im_per_spectrum_with_limits(
            bounds,
            ms_level,
            ExperimentMobilityLimits::default(),
        )
    }
    /// As [`Self::get_2d_peak_data_im_per_spectrum`], with explicit ceilings.
    pub fn get_2d_peak_data_im_per_spectrum_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        limits: ExperimentMobilityLimits,
    ) -> Result<SpectrumPeakDataIm> {
        let mut result = SpectrumPeakDataIm::default();
        self.append_2d_peak_data_im_per_spectrum_with_limits(
            bounds,
            ms_level,
            &mut result,
            limits,
        )?;
        Ok(result)
    }
    /// Append the selected rows to existing row vectors, keeping their contents.
    ///
    /// # Errors
    ///
    /// As [`Self::get_2d_peak_data_im_per_spectrum`]. Additionally returns
    /// [`Error::InvalidValue`] when the first selected retention time is
    /// exactly `-1` and the output holds no row to continue: the source's row
    /// cursor starts at `-1` (`MSExperiment.cpp:181`), so it would append
    /// through `mz.back()` into an empty vector.
    pub fn append_2d_peak_data_im_per_spectrum(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut SpectrumPeakDataIm,
    ) -> Result<()> {
        self.append_2d_peak_data_im_per_spectrum_with_limits(
            bounds,
            ms_level,
            output,
            ExperimentMobilityLimits::default(),
        )
    }
    /// As [`Self::append_2d_peak_data_im_per_spectrum`], with explicit ceilings.
    pub fn append_2d_peak_data_im_per_spectrum_with_limits(
        &self,
        bounds: AreaBounds,
        ms_level: usize,
        output: &mut SpectrumPeakDataIm,
        limits: ExperimentMobilityLimits,
    ) -> Result<()> {
        let mut work = Work::new(limits);
        let points = selected(self, bounds, ms_level, limits, &mut work)?;
        if points.len() == 0 {
            return Ok(());
        }
        let old_rows = output.rt.len();
        if old_rows != output.mz.len()
            || old_rows != output.intensity.len()
            || old_rows != output.ion_mobility.len()
        {
            return Err(invalid("bulk peak row vectors are not aligned"));
        }
        cap(old_rows, limits.max_output_rows)?;
        work.consume(old_rows)?;
        let mut old_points = 0usize;
        for ((mz, intensity), mobility) in output
            .mz
            .iter()
            .zip(&output.intensity)
            .zip(&output.ion_mobility)
        {
            if mz.len() != intensity.len() || mz.len() != mobility.len() {
                return Err(invalid("bulk peak row values are not aligned"));
            }
            old_points = add(old_points, mz.len())?;
            cap(old_points, limits.max_output_points)?;
        }
        // First pass: validate every conversion and mobility lookup and count
        // the rows the selection produces, before anything is allocated or
        // written. It tracks the row cursor exactly as the writing pass does,
        // so a lookup that would fail there fails here instead.
        work.consume(mul(points.len(), 5)?)?;
        let mut new_rows = 0usize;
        let mut previous_rt = ROW_RT_START;
        let mut row_mobility: Option<&[f32]> = None;
        for point in points.clone() {
            narrow(point.peak.mz, "m/z")?;
            finite_intensity(point.peak.intensity)?;
            let rt = narrow(point.spectrum.rt, "RT")?;
            if point.spectrum.rt != f64::from(previous_rt) {
                previous_rt = rt;
                row_mobility = frame_values(point.spectrum);
                new_rows = add(new_rows, 1)?;
                cap(add(old_rows, new_rows)?, limits.max_output_rows)?;
            } else if old_rows == 0 && new_rows == 0 {
                return Err(invalid(
                    "first selected RT is -1 but no existing output row exists",
                ));
            }
            if let Some(values) = row_mobility {
                value_at(values, point.peak_index)?;
            }
        }
        let rows = add(old_rows, new_rows)?;
        let points_total = add(old_points, points.len())?;
        cap(points_total, limits.max_output_points)?;
        // Old and new storage are charged together; only numeric payloads are
        // staged, and the caller's output is replaced only once it is complete.
        work.allocate::<f32>(rows)?;
        work.allocate::<Vec<f32>>(mul(rows, 3)?)?;
        work.allocate::<f32>(mul(points_total, 3)?)?;
        work.consume(add(mul(old_points, 3)?, mul(rows, 4)?)?)?;
        let mut staged = SpectrumPeakDataIm {
            rt: copied(&output.rt, rows, &mut work)?,
            mz: output.mz.clone(),
            intensity: output.intensity.clone(),
            ion_mobility: output.ion_mobility.clone(),
        };
        previous_rt = ROW_RT_START;
        row_mobility = None;
        for point in points {
            if point.spectrum.rt != f64::from(previous_rt) {
                previous_rt = point.spectrum.rt as f32;
                staged.rt.push(previous_rt);
                staged.mz.push(Vec::new());
                staged.intensity.push(Vec::new());
                staged.ion_mobility.push(Vec::new());
                // Source refreshes the IM array only when a row starts.
                row_mobility = frame_values(point.spectrum);
            }
            let mobility = match row_mobility {
                None => MOBILITY_NOT_SET,
                Some(values) => value_at(values, point.peak_index)?,
            };
            push_row(&mut staged.mz, point.peak.mz as f32)?;
            push_row(&mut staged.intensity, point.peak.intensity)?;
            push_row(&mut staged.ion_mobility, mobility)?;
        }
        *output = staged;
        Ok(())
    }

    /// Rasterize the run into a 2D intensity image (m/z by retention time).
    ///
    /// Bins peak intensities into a regular grid, producing the heatmap a
    /// viewer draws. The returned buffer holds `mz_bins * rt_bins` values in
    /// row-major order — m/z varies slowest, retention time fastest — so index
    /// `mz_bin * rt_bins + rt_bin` is one pixel. Rows are m/z bins, the y axis
    /// of a visualization; columns are retention-time bins, the x axis.
    ///
    /// Ports `rasterizeRTMZ` (`MSExperiment.cpp:262-518`). Only spectra of
    /// `RtMzRaster::ms_level` whose retention time lies in `[min_rt, max_rt]`
    /// contribute, and within each of them only the peaks in `[min_mz, max_mz]`.
    /// A peak exactly at `max_rt` or `max_mz` lands in the last bin rather than
    /// past the end, and intensities that share a pixel combine per the raster's
    /// `aggregation`.
    ///
    /// The source writes into a caller-allocated `float*` of exactly
    /// `rt_bins * mz_bins` entries and zero-fills it first; an owned [`Vec`]
    /// replaces both the pointer and the `Exception::NullPointer` the source
    /// throws for a null one. The source's `@note` recommending
    /// `numpy.empty` and its Python example describe that buffer and have no
    /// counterpart here. The source's other `@note`, that the run "should be
    /// sorted by RT and m/z", is a checked precondition in this port.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `rt_bins` or `mz_bins` is zero,
    /// when the pixel count overflows or exceeds [`Self::MAX_RASTER_PIXELS`],
    /// when the visited peaks exceed [`Self::MAX_MOBILITY_ITEMS`], or when a
    /// bound, retention time, m/z or intensity is not finite;
    /// [`Error::InvalidRange`] when `min_rt >= max_rt` or `min_mz >= max_mz`,
    /// as the source's two `Exception::InvalidRange` throws; and
    /// [`Error::UnsortedData`] when the spectra are not sorted by retention
    /// time or a contributing spectrum is not sorted by m/z.
    ///
    /// The source validates the bounds, then zero-fills the caller's buffer, and
    /// returns early for an empty run; here nothing is allocated before every
    /// check has passed. The source has no finiteness guard either: a NaN m/z
    /// passes both range tests and `static_cast<Int64>(NaN)` is undefined
    /// behaviour in C++.
    ///
    /// The source also carries two *negative-bin* guards this port has no
    /// counterpart for: it skips the whole spectrum when `rt_bin < 0`
    /// (`MSExperiment.cpp:362` in the single-threaded branch, `:425-428` in the
    /// parallel one) and skips one peak when `mz_bin < 0` (the `mz_bin >= 0`
    /// tests at `MSExperiment.cpp:376` and `:388`, and `:453`/`:472` in the
    /// parallel branch). Both are unreachable on the input the function
    /// documents: `RTBegin(min_rt)` and `MZBegin(min_mz)` are `lower_bound`
    /// calls, so every visited retention time is at least `min_rt` and every
    /// visited m/z at least `min_mz`, and the two scales are positive because
    /// `min < max` is enforced. They can only fire when the run is *not*
    /// sorted, which makes those `lower_bound` results meaningless — and that
    /// input this port rejects outright with [`Error::UnsortedData`] before any
    /// binning, through the `rt_begin`/`mz_begin` calls it shares with every
    /// other search. Nothing here therefore silently drops a spectrum or a
    /// peak; were the branch reachable, Rust's saturating `as usize` cast would
    /// clamp a negative bin into bin `0` rather than skip it.
    pub fn rasterize_rt_mz(&self, raster: &RtMzRaster) -> Result<Vec<f32>> {
        let pixels = raster.validate()?;
        self.cap_spectra()?;
        let begin = self.rt_begin(raster.min_rt)?;
        let end = self.rt_end(raster.max_rt)?;
        let mut windows = Vec::new();
        let mut visited = 0usize;
        for (index, spectrum) in self.spectra[begin..end].iter().enumerate() {
            if spectrum.ms_level != raster.ms_level {
                continue;
            }
            let first = spectrum.mz_begin(raster.min_mz)?;
            let last = spectrum.mz_end(raster.max_mz)?;
            visited = add(visited, last - first)?;
            if visited > Self::MAX_MOBILITY_ITEMS {
                return Err(invalid("visited peaks exceed MAX_MOBILITY_ITEMS"));
            }
            if first != last {
                windows
                    .try_reserve(1)
                    .map_err(|_| invalid("raster plan allocation failed"))?;
                windows.push((begin + index, first, last));
            }
        }
        for &(index, first, last) in &windows {
            let spectrum = &self.spectra[index];
            finite(spectrum.rt, "spectrum retention time")?;
            for peak in &spectrum.peaks[first..last] {
                finite(peak.mz, "peak m/z")?;
                finite(f64::from(peak.intensity), "peak intensity")?;
            }
        }
        let mut output = vec![0.0_f32; pixels];
        let rt_scale = raster.rt_bins as f64 / (raster.max_rt - raster.min_rt);
        let mz_scale = raster.mz_bins as f64 / (raster.max_mz - raster.min_mz);
        for &(index, first, last) in &windows {
            let spectrum = &self.spectra[index];
            let rt_bin = (((spectrum.rt - raster.min_rt) * rt_scale) as usize)
                .min(raster.rt_bins.saturating_sub(1));
            for peak in &spectrum.peaks[first..last] {
                let mz_bin = (((peak.mz - raster.min_mz) * mz_scale) as usize)
                    .min(raster.mz_bins.saturating_sub(1));
                let pixel = mz_bin * raster.rt_bins + rt_bin;
                match raster.aggregation {
                    RasterAggregation::Sum => output[pixel] += peak.intensity,
                    RasterAggregation::Max => {
                        if peak.intensity > output[pixel] {
                            output[pixel] = peak.intensity;
                        }
                    }
                }
            }
        }
        Ok(output)
    }

    fn cap_spectra(&self) -> Result<()> {
        if self.spectra.len() > Self::MAX_MOBILITY_ITEMS {
            return Err(invalid("spectrum count exceeds MAX_MOBILITY_ITEMS"));
        }
        Ok(())
    }
    fn check_mobility_order(&self) -> Result<()> {
        self.cap_spectra()?;
        check_sorted(&self.spectra, |spectrum| spectrum.drift_time)
    }
}

/// The per-peak ion mobility values of a spectrum, or `None` when it does not
/// present usable ones.
///
/// Source `maybeGetIMData` followed by the `unit != DriftTimeUnit::NONE` test
/// (`MSExperiment.cpp:194`, `MSExperiment.cpp:250`). The array-name rule itself
/// is [`MSSpectrum::maybe_im_data`], not re-derived here.
fn frame_values(spectrum: &MSSpectrum) -> Option<&[f32]> {
    match spectrum.maybe_im_data() {
        Some((unit, values)) if unit != DriftTimeUnit::None => Some(values),
        _ => None,
    }
}

/// Ion mobility of one peak of a spectrum, addressed by its index **within
/// that spectrum** as the source's `getPeakIndex().peak` is.
///
/// The source indexes the array without a bounds check, so an ion-mobility
/// array shorter than the peak list reads out of bounds; that is refused here.
fn peak_mobility(spectrum: &MSSpectrum, peak_index: usize) -> Result<f32> {
    match frame_values(spectrum) {
        None => Ok(MOBILITY_NOT_SET),
        Some(values) => value_at(values, peak_index),
    }
}
fn value_at(values: &[f32], index: usize) -> Result<f32> {
    let value = *values
        .get(index)
        .ok_or_else(|| invalid("ion mobility array is shorter than the exported peaks"))?;
    finite_intensity(value)
}

fn selected<'a>(
    experiment: &'a MSExperiment,
    bounds: AreaBounds,
    ms_level: usize,
    limits: ExperimentMobilityLimits,
    work: &mut Work,
) -> Result<AreaIter<'a>> {
    // Source Size -> UInt -> uint8_t -> int8_t -> UInt, with explicit modular
    // narrowing on all native platforms. Zero is not a wildcard.
    let options = AreaOptions::source_compatible(bounds, ms_level as u32);
    experiment.area_iter_with_budget(
        options,
        limits.max_spectra,
        limits.max_peaks,
        &mut work.remaining,
        &mut work.bytes,
    )
}
fn push_row(rows: &mut [Vec<f32>], value: f32) -> Result<()> {
    let row = rows
        .last_mut()
        .ok_or_else(|| invalid("bulk peak output has no current row"))?;
    row.try_reserve(1)
        .map_err(|_| invalid("bulk peak row allocation failed"))?;
    row.push(value);
    Ok(())
}
fn narrow(value: f64, label: &str) -> Result<f32> {
    let result = value as f32;
    if !value.is_finite() || !result.is_finite() {
        return Err(invalid(&format!(
            "bulk peak {label} cannot be represented as finite f32"
        )));
    }
    Ok(result)
}
fn finite_intensity(value: f32) -> Result<f32> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("bulk peak value is nonfinite"))
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    invalid("experiment mobility resource limit exceeded")
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
fn cap(count: usize, max: usize) -> Result<()> {
    if count > max { Err(limit()) } else { Ok(()) }
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Work {
    fn new(limits: ExperimentMobilityLimits) -> Self {
        Self {
            remaining: limits.max_work,
            bytes: limits.max_bytes,
        }
    }
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(count).ok_or_else(limit)?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<()> {
        if count != 0 {
            self.bytes = self
                .bytes
                .checked_sub(add(mul(count, size_of::<T>())?, 64)?)
                .ok_or_else(limit)?;
        }
        Ok(())
    }
}
fn copied<T: Copy>(input: &[T], capacity: usize, work: &mut Work) -> Result<Vec<T>> {
    work.allocate::<T>(capacity)?;
    let mut result = Vec::new();
    result.try_reserve_exact(capacity).map_err(|_| limit())?;
    result.extend_from_slice(input);
    Ok(result)
}
