// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ion mobility on a spectrum: drift time, IM data arrays, IM sorting and frame
//! rasterization.
//!
//! Source: the ion-mobility half of `KERNEL/MSSpectrum.h` /
//! `KERNEL/MSSpectrum.cpp` at Core SDK `bc9cc12`, together with the array-name
//! interpretation in `IONMOBILITY/IMDataArrayUtils.cpp` that
//! `MSSpectrum::containsIMData` delegates to. The support document
//! `docs/SPECTRUM_MOBILITY_SUPPORT.md` lists every source member, its
//! counterpart here, the preserved source conventions and the native
//! differences.
//!
//! Two representations coexist in the source and are both kept here. A
//! conventional spectrum carries one drift time for the whole scan, stored in
//! the `drift_time` field with `-1` meaning unset; an *IM frame* annotates every
//! peak, through a float data array whose *name* is a child of
//! `MS:1002893 ! ion mobility array` or one of three vendor `UserParam`
//! prefixes. [`MSSpectrum::contains_im_data`](crate::kernel::MSSpectrum::contains_im_data)
//! answers which one a spectrum is.
//!
//! Every operation whose cost scales with the peak count or the raster size is
//! checked against [`MSSpectrum::MAX_MOBILITY_ITEMS`](crate::kernel::MSSpectrum::MAX_MOBILITY_ITEMS)
//! or [`MSSpectrum::MAX_RASTER_PIXELS`](crate::kernel::MSSpectrum::MAX_RASTER_PIXELS)
//! before anything is allocated or mutated, so a rejected call leaves the
//! spectrum unchanged.
//!
//! The source parallelises neither of the sorts nor the rasterizer, so this
//! module has no OpenMP gap to record; the crate is serial throughout.

use super::{MSSpectrum, finite};
use crate::concept::constants::user_param;
use crate::error::{Error, Result};
use crate::metadata::DriftTimeUnit;

/// The sentinel the source uses for an unset drift time
/// (`IMTypes::DRIFTTIME_NOT_SET`, `IMTypes.h:93`).
const DRIFTTIME_NOT_SET: f64 = -1.0;

/// The nine PSI-MS children of `MS:1002893 ! ion mobility array` and the drift
/// time unit each one carries, pinned at the CV shipped in `resources/cv`.
///
/// The source resolves this at run time: `IMDataArrayUtils::getIMUnit` looks the
/// array name up in the PSI-MS controlled vocabulary, asks whether the term is a
/// child of `MS:1002893`, then reads the term's `has_units` relationship and maps
/// `MS:1002814` to [`DriftTimeUnit::InverseReducedMobility`], `UO:0000028` to
/// [`DriftTimeUnit::Millisecond`] and `UO:0000324` to
/// [`DriftTimeUnit::CollisionCrossSection`], in that order. At the pinned CV the
/// nine children resolve to the three inverse-reduced arrays (`MS:1003006`,
/// `MS:1003008`, `MS:1003155`) and six millisecond arrays; no child declares
/// `UO:0000324`, so the collision-cross-section branch is unreachable from the
/// CV and only the `UserParam` fallback below can produce it.
///
/// Pinning the table keeps this module free of the controlled-vocabulary
/// feature and of a CV file at run time. The name list is the same one
/// `crate::kernel::ranges` uses privately; the two must stay in step, and the
/// support document records that constraint.
const CV_IM_ARRAYS: [(&str, DriftTimeUnit); 9] = [
    (
        "mean ion mobility drift time array",
        DriftTimeUnit::Millisecond,
    ),
    ("mean ion mobility array", DriftTimeUnit::Millisecond),
    (
        "mean inverse reduced ion mobility array",
        DriftTimeUnit::InverseReducedMobility,
    ),
    ("raw ion mobility array", DriftTimeUnit::Millisecond),
    (
        "raw inverse reduced ion mobility array",
        DriftTimeUnit::InverseReducedMobility,
    ),
    (
        "raw ion mobility drift time array",
        DriftTimeUnit::Millisecond,
    ),
    (
        "deconvoluted ion mobility array",
        DriftTimeUnit::Millisecond,
    ),
    (
        "deconvoluted inverse reduced ion mobility array",
        DriftTimeUnit::InverseReducedMobility,
    ),
    (
        "deconvoluted ion mobility drift time array",
        DriftTimeUnit::Millisecond,
    ),
];

/// The unit an ion-mobility float data array name denotes, or `None` when the
/// name does not describe ion mobility at all.
///
/// This is `IMDataArrayUtils::getIMUnit` (`IMDataArrayUtils.cpp:37-93`): the
/// exact PSI-MS term first, then the two inverse-reduced `UserParam` prefixes,
/// then the generic `"Ion Mobility"` prefix whose embedded accession decides
/// between `1/K0`, CCS and milliseconds.
fn array_unit(name: &str) -> Option<DriftTimeUnit> {
    if let Some(&(_, unit)) = CV_IM_ARRAYS.iter().find(|(term, _)| *term == name) {
        return Some(unit);
    }
    if name.starts_with(user_param::MEAN_INVERSE_REDUCED_ION_MOBILITY_ARRAY)
        || name.starts_with(user_param::INVERSE_REDUCED_ION_MOBILITY)
    {
        return Some(DriftTimeUnit::InverseReducedMobility);
    }
    if name.starts_with(user_param::ION_MOBILITY) {
        if name.contains("MS:1002815") || name.contains("MS:1003006") {
            return Some(DriftTimeUnit::InverseReducedMobility);
        }
        if name.contains("MS:1002954") {
            return Some(DriftTimeUnit::CollisionCrossSection);
        }
        return Some(DriftTimeUnit::Millisecond);
    }
    None
}

/// First float data array that describes ion mobility, with its unit.
///
/// The source scans the arrays in storage order and stops at the first match
/// (`MSSpectrum.cpp:808-818`), so a spectrum carrying two ion-mobility arrays
/// uses the earlier one.
fn im_array(spectrum: &MSSpectrum) -> Option<(usize, DriftTimeUnit)> {
    spectrum
        .float_data_arrays
        .iter()
        .enumerate()
        .find_map(|(index, array)| array_unit(&array.name).map(|unit| (index, unit)))
}

/// How overlapping peaks combine inside one raster pixel.
///
/// Source `MSSpectrum::RasterAggregation` (`MSSpectrum.h:671-675`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RasterAggregation {
    /// Sum the intensities of all peaks falling into a pixel. The source
    /// default.
    #[default]
    Sum,
    /// Take the maximum intensity of all peaks falling into a pixel.
    Max,
}

/// The regular grid an ion-mobility frame is rasterized onto.
///
/// Source `MSSpectrum::rasterizeIMFrame`'s seven scalar parameters
/// (`MSSpectrum.h:720-728`), gathered into one value so the call site names each
/// bound. The m/z axis is the slow (row) axis and the ion-mobility axis the fast
/// (column) axis, matching the source's row-major output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImFrameRaster {
    /// Number of bins along the ion-mobility axis (image width). Must be
    /// positive.
    pub im_bins: usize,
    /// Number of bins along the m/z axis (image height). Must be positive.
    pub mz_bins: usize,
    /// Inclusive lower end of the ion-mobility range, in the array's own unit.
    pub min_im: f64,
    /// Inclusive upper end of the ion-mobility range; must exceed `min_im`.
    pub max_im: f64,
    /// Inclusive lower end of the m/z range, in Th.
    pub min_mz: f64,
    /// Inclusive upper end of the m/z range; must exceed `min_mz`.
    pub max_mz: f64,
    /// How peaks sharing a pixel combine.
    pub aggregation: RasterAggregation,
}

impl ImFrameRaster {
    /// A grid over the given ranges, summing peaks that share a pixel as the
    /// source default does.
    pub fn new(
        im_bins: usize,
        mz_bins: usize,
        min_im: f64,
        max_im: f64,
        min_mz: f64,
        max_mz: f64,
    ) -> Self {
        Self {
            im_bins,
            mz_bins,
            min_im,
            max_im,
            min_mz,
            max_mz,
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
    /// Returns [`Error::InvalidValue`] when `im_bins * mz_bins` overflows
    /// `usize`. The source multiplies unchecked (`MSSpectrum.cpp:895`) and then
    /// writes `im_bins * mz_bins` floats through a caller-supplied pointer, so
    /// an overflowing product there silently under-allocates.
    pub fn pixels(&self) -> Result<usize> {
        self.im_bins
            .checked_mul(self.mz_bins)
            .ok_or_else(|| Error::InvalidValue("raster pixel count overflows".into()))
    }

    fn validate(&self) -> Result<usize> {
        if self.im_bins == 0 {
            return Err(Error::InvalidValue(
                "number of IM bins must be positive".into(),
            ));
        }
        if self.mz_bins == 0 {
            return Err(Error::InvalidValue(
                "number of m/z bins must be positive".into(),
            ));
        }
        finite(self.min_im, "minimum ion mobility")?;
        finite(self.max_im, "maximum ion mobility")?;
        finite(self.min_mz, "minimum m/z")?;
        finite(self.max_mz, "maximum m/z")?;
        if self.min_im >= self.max_im {
            return Err(Error::InvalidRange(
                "minimum ion mobility must be below the maximum".into(),
            ));
        }
        if self.min_mz >= self.max_mz {
            return Err(Error::InvalidRange(
                "minimum m/z must be below the maximum".into(),
            ));
        }
        let pixels = self.pixels()?;
        if pixels > MSSpectrum::MAX_RASTER_PIXELS {
            return Err(Error::InvalidValue(
                "raster pixel count exceeds MAX_RASTER_PIXELS".into(),
            ));
        }
        Ok(pixels)
    }
}

/// A half-open run of peaks that is already known to be sorted, or known not to
/// be.
///
/// Source `MSSpectrum::Chunk` (`MSSpectrum.h:58-65`). Build these with
/// [`Chunks`] rather than by hand whenever the peaks are appended in runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Chunk {
    /// First peak index in the run, inclusive.
    pub start: usize,
    /// One past the last peak index in the run.
    pub end: usize,
    /// Whether the peaks in `[start, end)` are already sorted by m/z.
    pub is_sorted: bool,
}

impl Chunk {
    /// A run covering `[start, end)`, declared sorted or unsorted.
    pub const fn new(start: usize, end: usize, is_sorted: bool) -> Self {
        Self {
            start,
            end,
            is_sorted,
        }
    }
}

/// Accumulates [`Chunk`]s while peaks are appended to a spectrum.
///
/// Source `MSSpectrum::Chunks` (`MSSpectrum.h:77-91`). Each [`Chunks::add`]
/// closes a run that ends at the spectrum's current length and starts the next
/// one there, so the runs are contiguous and cover `[0, len)` once the caller
/// closes the last one.
///
/// The source holds a `const MSSpectrum&` for the whole life of the builder and
/// reads `size()` inside `add()`, which means the spectrum is mutated while a
/// reference to it is held. That cannot be expressed in Rust, so the spectrum is
/// passed to [`Chunks::add`] instead; the recorded chunks are identical.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Chunks {
    chunks: Vec<Chunk>,
}

impl Chunks {
    /// A builder with no runs recorded yet.
    pub fn new() -> Self {
        Self::default()
    }

    /// Close a run that ends at `spectrum`'s current peak count.
    ///
    /// The run starts where the previous one ended, or at 0 for the first.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the spectrum has shrunk below the
    /// previous run's end, or when the run count would exceed
    /// [`MSSpectrum::MAX_MOBILITY_ITEMS`]. The source computes `start > end` in
    /// that case and later forms an invalid iterator pair from it
    /// (`MSSpectrum.h:82`).
    pub fn add(&mut self, spectrum: &MSSpectrum, is_sorted: bool) -> Result<()> {
        if self.chunks.len() >= MSSpectrum::MAX_MOBILITY_ITEMS {
            return Err(Error::InvalidValue(
                "chunk count exceeds MAX_MOBILITY_ITEMS".into(),
            ));
        }
        let start = self.chunks.last().map_or(0, |chunk| chunk.end);
        let end = spectrum.len();
        if end < start {
            return Err(Error::InvalidValue(
                "spectrum shrank below the previous chunk end".into(),
            ));
        }
        self.chunks.push(Chunk::new(start, end, is_sorted));
        Ok(())
    }

    /// The recorded runs, in order (source `getChunks()`).
    pub fn chunks(&self) -> &[Chunk] {
        &self.chunks
    }

    /// Number of recorded runs.
    pub fn len(&self) -> usize {
        self.chunks.len()
    }

    /// Whether no run has been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.chunks.is_empty()
    }
}

/// Stable bottom-up merge of adjacent sorted runs of `indices`.
///
/// `bounds` holds the run boundaries, `bounds[0] == 0` and
/// `bounds[bounds.len() - 1] == indices.len()`. Merging adjacent stable runs
/// pairwise yields the same permutation as the source's balanced
/// `std::inplace_merge` recursion (`MSSpectrum.cpp:424-436`), because a stable
/// merge of adjacent runs is associative in its result.
fn merge_runs(indices: &[usize], bounds: &[usize], key: impl Fn(usize) -> f64) -> Vec<usize> {
    let mut runs: Vec<Vec<usize>> = bounds
        .windows(2)
        .map(|pair| indices[pair[0]..pair[1]].to_vec())
        .collect();
    while runs.len() > 1 {
        let mut merged: Vec<Vec<usize>> = Vec::with_capacity(runs.len().div_ceil(2));
        let mut pairs = runs.into_iter();
        while let Some(left) = pairs.next() {
            match pairs.next() {
                None => merged.push(left),
                Some(right) => {
                    let mut out = Vec::with_capacity(left.len() + right.len());
                    let (mut i, mut j) = (0, 0);
                    while i < left.len() && j < right.len() {
                        // `<` on the right element keeps the merge stable: an
                        // equal pair takes the left (earlier) element first.
                        if key(right[j]) < key(left[i]) {
                            out.push(right[j]);
                            j += 1;
                        } else {
                            out.push(left[i]);
                            i += 1;
                        }
                    }
                    out.extend_from_slice(&left[i..]);
                    out.extend_from_slice(&right[j..]);
                    merged.push(out);
                }
            }
        }
        runs = merged;
    }
    runs.pop().unwrap_or_default()
}

impl MSSpectrum {
    /// Largest peak, chunk or ion-mobility value count the operations in this
    /// module will process. Matches `RangeManager::MAX_ITEMS`, so a spectrum
    /// whose ranges can be computed can also be sorted by ion mobility.
    pub const MAX_MOBILITY_ITEMS: usize = 100_000_000;

    /// Largest raster the frame rasterizer will allocate, in pixels. At `f32`
    /// per pixel this caps one image at 64 MiB.
    pub const MAX_RASTER_PIXELS: usize = 16_777_216;

    /// The spectrum-wide ion mobility drift time, or `None` when it is unset.
    ///
    /// Drift times may be stored directly as an attribute of the spectrum, if
    /// they relate to the spectrum as a whole. For an ion mobility spectrum the
    /// drift time of the spectrum is always set here, while the drift time
    /// attribute of the [`Precursor`](crate::kernel::Precursor) may often be
    /// unpopulated.
    ///
    /// Source `getDriftTime` returns the raw `double`, with
    /// `IMTypes::DRIFTTIME_NOT_SET` (`-1`) meaning unset; the field
    /// `MSSpectrum::drift_time` keeps that representation, and this accessor
    /// turns the sentinel into `None` so callers cannot feed `-1` into a
    /// calculation by accident. The unit is
    /// [`MSSpectrum::drift_time_unit`].
    pub fn drift_time_if_set(&self) -> Option<f64> {
        (self.drift_time != DRIFTTIME_NOT_SET).then_some(self.drift_time)
    }

    /// Whether a spectrum-wide drift time is set, i.e. the stored value is not
    /// the source's `-1` sentinel.
    pub fn has_drift_time(&self) -> bool {
        self.drift_time != DRIFTTIME_NOT_SET
    }

    /// Set or clear the spectrum-wide ion mobility drift time.
    ///
    /// `None` writes the source sentinel `-1`. Setting the drift time does not
    /// touch [`MSSpectrum::drift_time_unit`],
    /// exactly as source `setDriftTime` does not.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `drift_time` is `Some` and not
    /// finite. The source assigns whatever it is given, so a NaN drift time
    /// reaches the range manager and every later comparison; this rejects it at
    /// the assignment and leaves the spectrum unchanged.
    pub fn set_drift_time(&mut self, drift_time: Option<f64>) -> Result<()> {
        match drift_time {
            None => self.drift_time = DRIFTTIME_NOT_SET,
            Some(value) => {
                finite(value, "spectrum drift time")?;
                self.drift_time = value;
            }
        }
        Ok(())
    }

    /// The drift time unit as the source spells it: `"<NONE>"`, `"ms"`,
    /// `"1/K0"`, `"FAIMS_CV"` or `"CCS"`.
    ///
    /// Source `getDriftTimeUnitAsString` indexes `NamesOfDriftTimeUnit`
    /// (`MSSpectrum.cpp:638-641`, `IMTypes.cpp:23`), which is the same table
    /// [`DriftTimeUnit::name`] carries. The source indexes without a bounds
    /// check, so a `DriftTimeUnit` holding `SIZE_OF_DRIFTTIMEUNIT` reads past
    /// the array; the Rust enum has no such member and the match is total.
    pub fn drift_time_unit_as_string(&self) -> &'static str {
        self.drift_time_unit.name()
    }

    /// Whether any float data array name marks this spectrum as an ion mobility
    /// frame, i.e. names a child of `MS:1002893 ! ion mobility array` or one of
    /// the three vendor `UserParam` prefixes.
    ///
    /// Source `containsIMData` (`MSSpectrum.cpp:820-825`). Only *names* are
    /// inspected; the array may be empty, or a different length from the peak
    /// list, and this still reports `true`. The operations that read the values
    /// check the length themselves.
    ///
    /// `crate::kernel::ranges` applies the same rule privately when it decides
    /// whether a spectrum's mobility range comes from an array or from the
    /// scalar drift time; the two must agree, and the support document records
    /// that.
    pub fn contains_im_data(&self) -> bool {
        im_array(self).is_some()
    }

    /// Index of the ion mobility float data array and the unit it carries.
    ///
    /// Source `getIMData` (`MSSpectrum.cpp:827-842`). Only works for spectra
    /// which represent an IM frame; test that first with
    /// [`MSSpectrum::contains_im_data`], or use
    /// [`MSSpectrum::maybe_im_data`] which does not fail.
    ///
    /// The unit can still be [`DriftTimeUnit::None`] when the array names a CV
    /// term that declares no usable unit. The source logs a warning in that
    /// case and returns the pair anyway; this port returns it silently, because
    /// the caller can see the unit.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when no float array names ion
    /// mobility data, as the source's `Exception::MissingInformation`.
    pub fn im_data(&self) -> Result<(usize, DriftTimeUnit)> {
        im_array(self).ok_or_else(|| {
            Error::MissingInformation(format!(
                "cannot get ion mobility data: no float array with a matching name among {}",
                self.float_data_arrays.len()
            ))
        })
    }

    /// The ion mobility values and their unit, or `None` when this spectrum is
    /// not an IM frame.
    ///
    /// Source `maybeGetIMData` (`MSSpectrum.cpp:844-855`) returns
    /// `{DriftTimeUnit::NONE, {}}` when there is no ion mobility array, which
    /// collides with the two cases that legitimately produce an empty vector or
    /// a `NONE` unit: an IM array with no entries, and an IM array whose CV term
    /// declares no unit. Returning an [`Option`] of a borrowed slice separates
    /// all three, and copies nothing.
    pub fn maybe_im_data(&self) -> Option<(DriftTimeUnit, &[f32])> {
        im_array(self).map(|(index, unit)| (unit, self.float_data_arrays[index].data.as_slice()))
    }

    /// Whether the peaks are ordered by their associated ion mobility value.
    ///
    /// Source `isSortedByIM` (`MSSpectrum.cpp:506-512`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when the spectrum carries no ion
    /// mobility array, as the source does.
    ///
    /// Returns [`Error::InvalidValue`] in two cases the source answers instead.
    /// When the array's length differs from the peak count the source reports on
    /// the array alone, so an empty ion mobility array makes an arbitrarily
    /// unsorted spectrum report *sorted*; the peaks it is supposed to describe
    /// are never consulted. And when the array holds a non-finite value, the
    /// source's `std::is_sorted` uses `<`, under which a NaN neighbour compares
    /// unordered and the array is again reported sorted. The port refuses both
    /// rather than returning a meaningless `true`. That is the one place this
    /// differs from [`MSSpectrum::is_sorted`], which folds an invalid coordinate
    /// into `false`.
    pub fn is_sorted_by_im(&self) -> Result<bool> {
        let values = self.checked_im_values()?;
        for &value in values {
            finite(f64::from(value), "ion mobility value")?;
        }
        Ok(values.windows(2).all(|pair| pair[0] <= pair[1]))
    }

    /// The ion mobility values, checked against the peak count.
    fn checked_im_values(&self) -> Result<&[f32]> {
        let (index, _) = self.im_data()?;
        let values = self.float_data_arrays[index].data.as_slice();
        if values.len() != self.peaks.len() {
            return Err(Error::InvalidValue(format!(
                "ion mobility array '{}' has {} entries for {} peaks",
                self.float_data_arrays[index].name,
                values.len(),
                self.peaks.len()
            )));
        }
        if values.len() > Self::MAX_MOBILITY_ITEMS {
            return Err(Error::InvalidValue(
                "ion mobility value count exceeds MAX_MOBILITY_ITEMS".into(),
            ));
        }
        Ok(values)
    }

    /// Sort the peaks, and every aligned data array, by ion mobility.
    ///
    /// Requires a float data array which is a child of
    /// `MS:1002893 ! ion mobility array`; see [`MSSpectrum::im_data`]. The sort
    /// is stable, so peaks sharing an ion mobility value keep their relative
    /// order, and it is skipped entirely when the array is already ordered, as
    /// in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when
    /// [`MSSpectrum::contains_im_data`] is false, as the source's
    /// `Exception::MissingInformation`. Returns [`Error::InvalidValue`] when the
    /// ion mobility array's length differs from the peak count, when it holds a
    /// non-finite value, when a different non-empty data array's length differs
    /// from the peak count, or when the peak count exceeds
    /// [`MSSpectrum::MAX_MOBILITY_ITEMS`].
    ///
    /// The length check is stricter than the source, deliberately. The source
    /// tests `std::is_sorted` on the ion mobility array before it sorts
    /// (`MSSpectrum.cpp:390`), and an empty range is sorted, so a spectrum whose
    /// ion mobility array has no entries — which `checkDataArraySizes_`
    /// explicitly permits — is silently left in whatever order it was in, with
    /// no error and no indication that no peak carried a mobility value. This
    /// port requires one value per peak. Every check runs before any peak
    /// moves.
    pub fn sort_by_ion_mobility(&mut self) -> Result<()> {
        let values = self.checked_im_values()?;
        for &value in values {
            finite(f64::from(value), "ion mobility value")?;
        }
        if values.windows(2).all(|pair| pair[0] <= pair[1]) {
            return Ok(());
        }
        self.validate_data_arrays()?;
        let mut indices: Vec<usize> = (0..self.peaks.len()).collect();
        let (index, _) = self.im_data()?;
        let values = &self.float_data_arrays[index].data;
        indices.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).expect("finite"));
        self.select(&indices)
    }

    /// Sort the peaks by m/z, exploiting runs that are already sorted.
    ///
    /// `chunks` describes consecutive runs of the peak list; each says whether
    /// its peaks are already ordered. Unordered runs are sorted individually and
    /// all runs are then merged, which is cheaper than sorting the whole list
    /// when the input arrives in ordered batches. Build `chunks` with
    /// [`Chunks`]. The result is the stable m/z order, identical to
    /// [`MSSpectrum::sort_by_position`]; every aligned data array follows.
    ///
    /// An empty `chunks` sorts nothing and succeeds, as in the source
    /// (`MSSpectrum.cpp:399`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the runs do not tile `[0, len)`
    /// exactly — the first must start at 0, each must start where the previous
    /// ended, none may end before it starts, and the last must end at the peak
    /// count — when a non-empty data array's length differs from the peak count,
    /// when a peak m/z is not finite, or when the peak count exceeds
    /// [`MSSpectrum::MAX_MOBILITY_ITEMS`]. Returns [`Error::UnsortedData`] when
    /// a run marked `is_sorted` is not in fact sorted by m/z.
    ///
    /// Neither condition is checked by the source, and both are observable
    /// there. A run list that stops short of the peak count leaves the tail
    /// unsorted when the spectrum has data arrays but sorts the whole list when
    /// it has none, because the two branches treat `chunks` differently
    /// (`MSSpectrum.cpp:405-441`); and a run that lies about being sorted feeds
    /// `std::inplace_merge` an unsorted range, whose result is unspecified.
    /// Checking costs one linear scan, which the merge pays anyway.
    pub fn sort_by_position_presorted(&mut self, chunks: &[Chunk]) -> Result<()> {
        if chunks.is_empty() {
            return Ok(());
        }
        if self.peaks.len() > Self::MAX_MOBILITY_ITEMS {
            return Err(Error::InvalidValue(
                "peak count exceeds MAX_MOBILITY_ITEMS".into(),
            ));
        }
        if chunks.len() > Self::MAX_MOBILITY_ITEMS {
            return Err(Error::InvalidValue(
                "chunk count exceeds MAX_MOBILITY_ITEMS".into(),
            ));
        }
        for peak in &self.peaks {
            finite(peak.mz, "peak m/z")?;
        }
        self.validate_data_arrays()?;

        let mut bounds = Vec::with_capacity(chunks.len() + 1);
        bounds.push(0usize);
        for chunk in chunks {
            if chunk.start != *bounds.last().expect("seeded") {
                return Err(Error::InvalidValue(
                    "chunks must tile the peak list without gaps or overlaps".into(),
                ));
            }
            if chunk.end < chunk.start {
                return Err(Error::InvalidValue("chunk ends before it starts".into()));
            }
            bounds.push(chunk.end);
        }
        if *bounds.last().expect("seeded") != self.peaks.len() {
            return Err(Error::InvalidValue(
                "chunks must cover every peak of the spectrum".into(),
            ));
        }
        for chunk in chunks {
            if chunk.is_sorted
                && self.peaks[chunk.start..chunk.end]
                    .windows(2)
                    .any(|pair| pair[0].mz > pair[1].mz)
            {
                return Err(Error::UnsortedData);
            }
        }
        if chunks.len() == 1 && chunks[0].is_sorted {
            return Ok(());
        }

        let mut indices: Vec<usize> = (0..self.peaks.len()).collect();
        for chunk in chunks {
            if !chunk.is_sorted {
                indices[chunk.start..chunk.end].sort_by(|&a, &b| {
                    self.peaks[a]
                        .mz
                        .partial_cmp(&self.peaks[b].mz)
                        .expect("finite")
                });
            }
        }
        let merged = merge_runs(&indices, &bounds, |index| self.peaks[index].mz);
        self.select(&merged)
    }

    /// Rasterize an ion mobility frame into a 2D intensity image (m/z by ion
    /// mobility).
    ///
    /// Bins the peak intensities into a regular grid, producing a heatmap of an
    /// IM frame: a spectrum in which every peak carries an ion mobility value in
    /// a float data array. The returned buffer holds `mz_bins * im_bins` values
    /// in row-major order — m/z varies slowest, ion mobility fastest — so index
    /// `mz_bin * im_bins + im_bin` is one pixel. Rows are m/z bins, the y axis
    /// of a visualization; columns are ion mobility bins, the x axis.
    ///
    /// Peaks outside the requested ranges are skipped. A peak exactly at
    /// `max_mz` or `max_im` lands in the last bin rather than past the end, and
    /// intensities that share a pixel combine per
    /// the raster's `aggregation`.
    ///
    /// Unlike a retention-time raster this does not require the spectrum to be
    /// sorted by m/z, because IM frame data is typically unsorted; all peaks are
    /// visited linearly.
    ///
    /// The source writes into a caller-allocated `float*` of exactly
    /// `im_bins * mz_bins` entries and zero-initializes it first; an owned
    /// [`Vec`] replaces both the pointer and the `Exception::NullPointer` the
    /// source throws for a null one. The source's Python example is not
    /// reproduced here, because it describes the NumPy buffer the pyOpenMS
    /// binding passes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::MissingInformation`] when the spectrum has no ion
    /// mobility array; [`Error::InvalidValue`] when `im_bins` or `mz_bins` is
    /// zero, when the pixel count overflows or exceeds
    /// [`MSSpectrum::MAX_RASTER_PIXELS`], when the peak count exceeds
    /// [`MSSpectrum::MAX_MOBILITY_ITEMS`], when the ion mobility array's length
    /// differs from the peak count, or when a bound, peak or ion mobility value
    /// is not finite; and [`Error::InvalidRange`] when `min_im >= max_im` or
    /// `min_mz >= max_mz`.
    ///
    /// The source validates the bounds, then zeroes the caller's buffer, and
    /// only then compares the ion mobility array's length against the peak count
    /// (`MSSpectrum.cpp:891-906`) — so that particular rejection still clears
    /// the caller's image. Here nothing is allocated before every check has
    /// passed. The source also has no non-finite guard: a NaN m/z passes both
    /// range tests, and `static_cast<Int64>(NaN)` is undefined behaviour in C++.
    pub fn rasterize_im_frame(&self, raster: &ImFrameRaster) -> Result<Vec<f32>> {
        let pixels = raster.validate()?;
        if self.peaks.len() > Self::MAX_MOBILITY_ITEMS {
            return Err(Error::InvalidValue(
                "peak count exceeds MAX_MOBILITY_ITEMS".into(),
            ));
        }
        // Presence is required even for an empty spectrum, as in the source.
        self.im_data()?;
        if self.peaks.is_empty() {
            return Ok(vec![0.0; pixels]);
        }
        let values = self.checked_im_values()?;
        for peak in &self.peaks {
            finite(peak.mz, "peak m/z")?;
            finite(f64::from(peak.intensity), "peak intensity")?;
        }
        for &value in values {
            finite(f64::from(value), "ion mobility value")?;
        }

        let mut output = vec![0.0_f32; pixels];
        let im_scale = raster.im_bins as f64 / (raster.max_im - raster.min_im);
        let mz_scale = raster.mz_bins as f64 / (raster.max_mz - raster.min_mz);
        for (peak, &mobility) in self.peaks.iter().zip(values) {
            let mobility = f64::from(mobility);
            if peak.mz < raster.min_mz
                || peak.mz > raster.max_mz
                || mobility < raster.min_im
                || mobility > raster.max_im
            {
                continue;
            }
            let mz_bin = (((peak.mz - raster.min_mz) * mz_scale) as usize).min(raster.mz_bins - 1);
            let im_bin = (((mobility - raster.min_im) * im_scale) as usize).min(raster.im_bins - 1);
            let pixel = mz_bin * raster.im_bins + im_bin;
            match raster.aggregation {
                RasterAggregation::Sum => output[pixel] += peak.intensity,
                RasterAggregation::Max => {
                    if peak.intensity > output[pixel] {
                        output[pixel] = peak.intensity;
                    }
                }
            }
        }
        Ok(output)
    }
}
