// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Seeds, mass traces and isotope patterns of the centroided-peptide feature
//! finder (`FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h`).
//!
//! The source wraps five helper types in one stateless struct,
//! `FeatureFinderAlgorithmPickedHelperStructs`, for `FeatureFinderAlgorithmPicked`
//! and the trace fitters `TraceFitter`, `GaussTraceFitter` and `EGHTraceFitter`
//! (not ported yet). The wrapper becomes this module. The complete API mapping,
//! the native differences and the evidence are in
//! `docs/FEATURE_FINDER_PICKED_HELPER_STRUCTS_SUPPORT.md`.
//!
//! # Peaks are indices plus copied values
//!
//! The source `MassTrace` stores `(retention time, const Peak1D*)` pairs whose
//! pointers reach into the peak map that `FeatureFinderAlgorithmPicked` owns.
//! Holding references instead would borrow that map for as long as any trace
//! lives, so a [`TracePeak`](crate::analysis::feature_finder_picked::helper_structs::TracePeak)
//! records the peak's `(spectrum, peak)` indices together with copies of the
//! retention time, m/z and intensity that the algorithm and the fitters read. A
//! copy is a snapshot. It agrees with the source's read through the pointer as
//! long as the map does not change while traces exist:
//! `FeatureFinderAlgorithmPicked` sorts the map before it builds a trace and
//! afterwards changes only float data arrays, never a peak's m/z or intensity.
//! The pointer identity the source compares (`max_peak == &peak`) becomes
//! equality of the `(spectrum, peak)` indices.
//!
//! # Precision
//!
//! Intensities stay `f32`, because the source `Peak1D::IntensityType` is
//! `float`. `MassTrace::avg_mz`, `MassTraces::update_baseline` and
//! `MassTraces::intensity_profile` promote each `f32` to `f64` exactly where the
//! source converts a `float` into a `double` expression, so their sums
//! accumulate in `f64`, not in `f32`. All three promote with the Linux x86_64
//! Release build's `cvtss2sd` (`scoring::x86_64::widen`) rather than with a
//! Rust cast, and `avg_mz` and `intensity_profile` also follow the executed
//! operand order of the additions, multiplications and the division, so a NaN
//! any of them creates or passes on carries the executed sign and payload on
//! every host.
//!
//! # Bounded work
//!
//! Operations that allocate in proportion to their input check a ceiling before
//! allocating: `IsotopePattern::MAX_SIZE`, `MassTraces::MAX_TRACES`,
//! `MassTraces::MAX_PEAKS` and `MassTraces::MAX_PROFILE_STEPS`. Every other
//! operation reads the traces once and allocates nothing. The source is serial,
//! and so is this module.

use std::ops::{Index, IndexMut};

use crate::kernel::{ConvexHull2D, Point2D};
use crate::{Error, Result};

/// The refusal of [`MassTraces::intensity_profile`] where a NaN retention time
/// meets the merge that never ends in the source (`CPP-242`).
pub(crate) const NAN_RT_MERGE_WHAT: &str =
    "a NaN retention time cannot be merged into an intensity profile";

/// Seed of a feature: a local intensity maximum in one spectrum (source
/// `FeatureFinderAlgorithmPickedHelperStructs::Seed`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Seed {
    /// Spectrum index.
    pub spectrum: usize,
    /// Peak index within the spectrum.
    pub peak: usize,
    /// Intensity, `f32` like the source `float`.
    pub intensity: f32,
}

impl Seed {
    /// A seed at the given spectrum and peak index with the given intensity.
    pub const fn new(spectrum: usize, peak: usize, intensity: f32) -> Self {
        Self {
            spectrum,
            peak,
            intensity,
        }
    }

    /// Whether this seed is less intense than `rhs`: source `Seed::operator<`.
    ///
    /// Compares only the `f32` intensities, with `<`. Seeds of equal intensity
    /// are therefore unordered whatever their positions, a NaN intensity is
    /// neither less nor greater than any other, and `-0.0` is not less than
    /// `+0.0`.
    ///
    /// The source overloads `operator<` for two uses in
    /// `FeatureFinderAlgorithmPicked`. `std::sort` orders the seeds by
    /// descending intensity with it. In debug mode it also orders the
    /// `std::map<Seed, std::string>` of abort reasons, where seeds of equal
    /// intensity are one key: the first seed stored keeps its position and takes
    /// the last reason. This port does not implement [`PartialOrd`], because an
    /// ordering that ignores the positions would contradict the structural
    /// [`PartialEq`] that compares them. The sorting policy, including ties, and
    /// the keying of an abort-reason map belong to the caller.
    pub fn is_less_intense_than(&self, rhs: &Self) -> bool {
        self.intensity < rhs.intensity
    }
}

/// One peak of a [`MassTrace`]: where it came from and the values the algorithm
/// reads.
///
/// Stands in for the source pair `(double, const Peak1D*)`; see the module
/// documentation. The indices are provenance only. Nothing in this module
/// dereferences them, and a caller whose peaks do not come from an experiment,
/// such as a class test, may choose any values.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TracePeak {
    /// Index of the spectrum the peak belongs to.
    pub spectrum: usize,
    /// Index of the peak within that spectrum.
    pub peak: usize,
    /// Retention time of the spectrum in seconds; the source pair's `first`.
    pub rt: f64,
    /// m/z of the peak, which the source reads through the pointer.
    pub mz: f64,
    /// Intensity of the peak, `f32` like `Peak1D::IntensityType`.
    pub intensity: f32,
}

impl TracePeak {
    /// A peak at `(spectrum, peak)` with the given retention time, m/z and
    /// intensity.
    pub const fn new(spectrum: usize, peak: usize, rt: f64, mz: f64, intensity: f32) -> Self {
        Self {
            spectrum,
            peak,
            rt,
            mz,
            intensity,
        }
    }
}

/// A mass trace: the peaks of one isotope over consecutive spectra (source
/// `FeatureFinderAlgorithmPickedHelperStructs::MassTrace`).
///
/// This is not [`crate::kernel::MassTrace`], the port of `KERNEL/MassTrace.h`
/// that mass-trace detection uses; the two source types share only a name.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassTrace {
    /// The maximum peak, set by [`Self::update_maximum`]; `None` until then, as
    /// the source pointer defaults to `nullptr`.
    pub max_peak: Option<TracePeak>,
    /// Retention time of the maximum peak.
    ///
    /// The source leaves this `double` uninitialised until
    /// [`Self::update_maximum`] sets it; it starts at `0.0` here.
    pub max_rt: f64,
    /// Theoretical intensity contribution of the trace, scaled to `[0, 1]`.
    ///
    /// Uninitialised in the source until the algorithm assigns it; `0.0` here.
    pub theoretical_int: f64,
    /// Contained peaks, in the order they were added.
    pub peaks: Vec<TracePeak>,
}

impl MassTrace {
    /// The convex hull of the trace's `(retention time, m/z)` points: source
    /// `getConvexhull`.
    ///
    /// Peaks with equal retention times share one scan whose m/z interval covers
    /// them all, so the peak order does not matter. The result is OpenMS's
    /// per-scan hull, as [`ConvexHull2D`] documents, not a geometric convex hull.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the trace holds more than
    /// [`MassTraces::MAX_PEAKS`] peaks, checked before allocating, or when a
    /// retention time or m/z is not finite. The source adds such values to a
    /// `std::map` keyed by retention time without a check, and a NaN key breaks
    /// that map's ordering.
    pub fn convex_hull(&self) -> Result<ConvexHull2D> {
        preflight_peaks(self.peaks.len())?;
        let points: Vec<Point2D> = self
            .peaks
            .iter()
            .map(|peak| Point2D::new(peak.rt, peak.mz))
            .collect();
        ConvexHull2D::from_points(&points)
    }

    /// Sets [`Self::max_peak`] and [`Self::max_rt`] to the most intense peak:
    /// source `updateMaximum`.
    ///
    /// Compares the `f32` intensities with a strict `>`. The first of several
    /// equally intense peaks therefore wins, a NaN intensity never replaces the
    /// current maximum, and a leading NaN is never replaced. An empty trace
    /// leaves both fields unchanged.
    pub fn update_maximum(&mut self) {
        let Some((first, rest)) = self.peaks.split_first() else {
            return;
        };
        let mut max = *first;
        for peak in rest {
            if peak.intensity > max.intensity {
                max = *peak;
            }
        }
        self.max_rt = max.rt;
        self.max_peak = Some(max);
    }

    /// The intensity-weighted average m/z of the trace: source `getAvgMZ`.
    ///
    /// Sums `mz * intensity` and `intensity` in `f64`, in peak order, promoting
    /// each `f32` intensity first, and divides the two sums. The quotient is
    /// the plain IEEE result, as in the source. An empty trace, or one whose
    /// intensities are all zero, yields NaN, and a NaN average never matches a
    /// seed in [`MassTraces::is_valid`].
    ///
    /// Every operation follows the Linux x86_64 Release build's SSE
    /// instructions (`libOpenMS.so` `0x18f1570`: `cvtss2sd` of the intensity,
    /// `addsd` onto the intensity sum, `mulsd` of the m/z by the intensity,
    /// `addsd` onto the product sum, `divsd` of the two sums), so a NaN
    /// carries the executed sign and payload on every host: `0 / 0` is x86_64's
    /// default NaN, whose sign bit is set, and the `.plot` file of
    /// `writeFeatureDebugInfo_` prints it as `-nan` (executed:
    /// `../oracle/ffap-complete-fix5`, a trace of zero intensities).
    pub fn avg_mz(&self) -> f64 {
        use crate::analysis::feature_finder_picked::scoring::x86_64;
        let mut sum = 0.0;
        let mut intensities = 0.0;
        for peak in &self.peaks {
            let intensity = x86_64::widen(peak.intensity);
            intensities = x86_64::add(intensities, intensity);
            sum = x86_64::add(sum, x86_64::mul(peak.mz, intensity));
        }
        x86_64::div(sum, intensities)
    }

    /// Whether the trace holds at least three peaks: source `isValid`, whose
    /// documentation says "more than 2 points".
    pub fn is_valid(&self) -> bool {
        self.peaks.len() >= 3
    }
}

/// The mass traces of one feature candidate (source
/// `FeatureFinderAlgorithmPickedHelperStructs::MassTraces`).
///
/// The source inherits privately from `std::vector<MassTrace>` and re-exports
/// only `size`, `at`, `reserve`, `push_back`, `operator[]`, `back`, `clear`,
/// `begin` and `end`. They map to [`Self::len`], [`Self::get`],
/// [`Self::reserve`], [`Self::push`], [`Index`], [`Self::last`],
/// [`Self::clear`] and [`Self::iter`], plus their mutable counterparts. As in
/// the source, no other vector operation, such as erase, insert or sort, is
/// offered.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassTraces {
    traces: Vec<MassTrace>,
    /// Index of the maximum-intensity trace; `0` by default, as the source
    /// constructor sets.
    pub max_trace: usize,
    /// Estimated baseline intensity in the region of the feature, used by the
    /// fit.
    ///
    /// The source leaves this `double` uninitialised until
    /// [`Self::update_baseline`] or the algorithm assigns it; it starts at `0.0`
    /// here.
    pub baseline: f64,
}

impl MassTraces {
    /// Most traces [`Self::reserve`] pre-allocates room for.
    pub const MAX_TRACES: usize = 100_000;

    /// Most peaks, over all traces, that [`Self::intensity_profile`] and
    /// [`MassTrace::convex_hull`] allocate for.
    ///
    /// A feature candidate holds a few isotope traces over at most a few
    /// thousand spectra, several orders of magnitude below this.
    pub const MAX_PEAKS: usize = 1_000_000;

    /// Most merge steps [`Self::intensity_profile`] may take, bounded before it
    /// starts.
    ///
    /// Merging trace `t` walks at most the profile built from the traces before
    /// it plus trace `t` itself, so the bound is the sum of those lengths over
    /// every trace after the first.
    pub const MAX_PROFILE_STEPS: usize = 100_000_000;

    /// An empty collection with `max_trace` 0 and baseline 0: source
    /// `MassTraces()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// The number of traces: source `size`.
    pub fn len(&self) -> usize {
        self.traces.len()
    }

    /// Whether the collection holds no traces.
    pub fn is_empty(&self) -> bool {
        self.traces.is_empty()
    }

    /// The trace at `index`, or `None` when out of range: source `at`, which
    /// throws `std::out_of_range` instead.
    pub fn get(&self, index: usize) -> Option<&MassTrace> {
        self.traces.get(index)
    }

    /// The mutable trace at `index`, or `None` when out of range.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut MassTrace> {
        self.traces.get_mut(index)
    }

    /// Pre-allocates room for at least `capacity` traces in total: source
    /// `reserve`.
    ///
    /// As with `std::vector::reserve`, the argument is the total capacity, not
    /// an additional one, and a capacity at or below the current length changes
    /// nothing. The ceiling bounds only the requested pre-allocation; traces
    /// added with [`Self::push`] grow the collection like a vector.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `capacity` exceeds
    /// [`Self::MAX_TRACES`] or the allocation fails, leaving the collection
    /// unchanged. The source fails only at the allocator's limit, by throwing
    /// `std::length_error` or `std::bad_alloc`.
    pub fn reserve(&mut self, capacity: usize) -> Result<()> {
        if capacity > Self::MAX_TRACES {
            return Err(Error::InvalidValue(format!(
                "reserving {capacity} mass traces exceeds the limit of {}",
                Self::MAX_TRACES
            )));
        }
        let additional = capacity.saturating_sub(self.traces.len());
        self.traces.try_reserve(additional).map_err(|error| {
            Error::InvalidValue(format!("cannot reserve {capacity} mass traces: {error}"))
        })
    }

    /// Appends a trace: source `push_back`.
    pub fn push(&mut self, trace: MassTrace) {
        self.traces.push(trace);
    }

    /// The last trace, or `None` when empty: source `back`, which is undefined
    /// on an empty vector.
    pub fn last(&self) -> Option<&MassTrace> {
        self.traces.last()
    }

    /// The mutable last trace, or `None` when empty.
    pub fn last_mut(&mut self) -> Option<&mut MassTrace> {
        self.traces.last_mut()
    }

    /// Removes every trace: source `clear`.
    ///
    /// [`Self::max_trace`] and [`Self::baseline`] keep their values, because the
    /// source's `std::vector::clear` does not reach them.
    pub fn clear(&mut self) {
        self.traces.clear();
    }

    /// The traces in order: source `begin`/`end`.
    pub fn iter(&self) -> std::slice::Iter<'_, MassTrace> {
        self.traces.iter()
    }

    /// The traces in order, mutably.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, MassTrace> {
        self.traces.iter_mut()
    }

    /// The traces as a slice.
    pub fn as_slice(&self) -> &[MassTrace] {
        &self.traces
    }

    /// The number of peaks over all traces: source `getPeakCount`.
    pub fn peak_count(&self) -> usize {
        self.traces.iter().map(|trace| trace.peaks.len()).sum()
    }

    /// Whether the collection still describes the seed: source
    /// `isValid(seed_mz, trace_tolerance)`.
    ///
    /// True when there are at least two traces and the [`MassTrace::avg_mz`] of
    /// some trace lies within `trace_tolerance` of `seed_mz`, inclusive:
    /// `|seed_mz - avg_mz| <= trace_tolerance`. The source method is not
    /// `const` but changes nothing, so this takes `&self`.
    ///
    /// # Arguments
    ///
    /// * `seed_mz` - m/z of the seed the traces were extended from.
    /// * `trace_tolerance` - largest accepted m/z distance between the seed and
    ///   a trace's average m/z.
    pub fn is_valid(&self, seed_mz: f64, trace_tolerance: f64) -> bool {
        if self.traces.len() < 2 {
            return false;
        }
        self.traces
            .iter()
            .any(|trace| (seed_mz - trace.avg_mz()).abs() <= trace_tolerance)
    }

    /// The index of the trace with the highest theoretical intensity: source
    /// `getTheoreticalmaxPosition`.
    ///
    /// Compares [`MassTrace::theoretical_int`] with a strict `>`, so the first
    /// of several equal maxima wins and a leading NaN is never replaced.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there are no traces, where the
    /// source throws `Exception::Precondition` in every build mode.
    pub fn theoretical_max_position(&self) -> Result<usize> {
        let Some((first, rest)) = self.traces.split_first() else {
            return Err(precondition(
                "There must be at least one trace to determine the theoretical maximum trace!",
            ));
        };
        let mut max = 0;
        let mut max_int = first.theoretical_int;
        for (offset, trace) in rest.iter().enumerate() {
            if trace.theoretical_int > max_int {
                max_int = trace.theoretical_int;
                max = offset + 1;
            }
        }
        Ok(max)
    }

    /// Sets [`Self::baseline`] to the lowest peak intensity over all traces:
    /// source `updateBaseline`.
    ///
    /// Without traces the baseline becomes `0.0`. Otherwise the first peak in
    /// trace order sets it and every later peak with a strictly lower intensity
    /// replaces it, compared in `f64` after promoting the `f32` with the
    /// Release build's `cvtss2sd` (`x86_64::widen`), so a NaN intensity gives
    /// the baseline the executed sign and payload on every host. A NaN first
    /// peak therefore leaves a NaN baseline, and a later NaN is skipped. When
    /// traces exist but none holds a peak, the baseline keeps its value, as in
    /// the source, where that value may still be uninitialised.
    ///
    /// The baseline is not debug-only: `run_` scales it (`.cpp:661`) and both
    /// fitters add it to every theoretical intensity (`.cpp:1983`, `.cpp:2098`)
    /// before the crop and quality correlations, so its bits reach the stored
    /// `score_fit` and `score_correlation` and the `.plot` formula.
    pub fn update_baseline(&mut self) {
        use crate::analysis::feature_finder_picked::scoring::x86_64;
        if self.traces.is_empty() {
            self.baseline = 0.0;
            return;
        }
        let mut first = true;
        for peak in self.traces.iter().flat_map(|trace| trace.peaks.iter()) {
            let intensity = x86_64::widen(peak.intensity);
            if first {
                self.baseline = intensity;
                first = false;
            }
            if intensity < self.baseline {
                self.baseline = intensity;
            }
        }
    }

    /// The smallest and largest peak retention time over all traces: source
    /// `getRTBounds`.
    ///
    /// Starts from `(f64::MAX, -f64::MAX)` and narrows with strict comparisons,
    /// so NaN retention times are skipped and traces without peaks leave those
    /// start values, as in the source.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there are no traces, where the
    /// source throws `Exception::Precondition` in every build mode.
    pub fn rt_bounds(&self) -> Result<(f64, f64)> {
        if self.traces.is_empty() {
            return Err(precondition(
                "There must be at least one trace to determine the RT boundaries!",
            ));
        }
        let mut min = f64::MAX;
        let mut max = -f64::MAX;
        for peak in self.traces.iter().flat_map(|trace| trace.peaks.iter()) {
            if peak.rt > max {
                max = peak.rt;
            }
            if peak.rt < min {
                min = peak.rt;
            }
        }
        Ok((min, max))
    }

    /// A flat representation of the traces: one summed intensity per retention
    /// time, comparable to the traces' total ion chromatogram. Source
    /// `computeIntensityProfile`.
    ///
    /// Each entry is `(retention time, intensity)`. The first trace is copied
    /// unchanged, and each later trace is merged into the profile in one
    /// forward walk. At the current profile entry, a trace peak with a smaller
    /// retention time is inserted before the entry; a peak with an exactly equal
    /// retention time is added to the entry, in `f64` after promoting its
    /// `f32` intensity, and both advance; a peak with a larger retention time
    /// moves the walk past the entry; peaks beyond the profile's end are
    /// appended. Sorted traces thus give a sorted merge. Unsorted traces, or
    /// repeated retention times within one trace, give the same possibly
    /// unsorted profile as the source's `std::list` walk.
    ///
    /// The source fills an output `std::list` that its documentation requires
    /// to be empty; returning a new vector is that contract. The list becomes an
    /// index-linked list internally, so each insertion costs what it costs in
    /// the source rather than shifting a vector.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`], before allocating, when the traces hold
    /// more than [`Self::MAX_PEAKS`] peaks or the merge could take more than
    /// [`Self::MAX_PROFILE_STEPS`] steps.
    ///
    /// Also returns [`Error::InvalidValue`] when a NaN retention time meets a
    /// profile entry during a merge. Every comparison is false there, so the
    /// source loop neither advances nor consumes a peak and never ends. A NaN
    /// that is only copied from the first trace or appended passes through, as
    /// in the source.
    ///
    /// An empty collection yields an empty profile. The source dereferences
    /// the first trace without a check, which is undefined behaviour when there
    /// is none.
    pub fn intensity_profile(&self) -> Result<Vec<(f64, f64)>> {
        let peaks = self.preflight_profile()?;
        let Some((first, rest)) = self.traces.split_first() else {
            return Ok(Vec::new());
        };
        // `cvtss2sd` of each intensity and, where two traces meet, `addsd` with
        // the new intensity as the destination (`libOpenMS.so` `0x18f1912`),
        // so a NaN (`inf - inf`, or two NaN operands) carries the Release
        // build's bits on every host.
        use crate::analysis::feature_finder_picked::scoring::x86_64;
        let mut profile = LinkedProfile::with_capacity(peaks);
        let mut previous = NIL;
        for peak in &first.peaks {
            previous = profile.insert_after(previous, (peak.rt, x86_64::widen(peak.intensity)));
        }
        for trace in rest {
            let mut previous = NIL;
            let mut current = profile.head;
            let mut index = 0;
            while let Some(peak) = trace.peaks.get(index) {
                let intensity = x86_64::widen(peak.intensity);
                if current == NIL {
                    previous = profile.insert_after(previous, (peak.rt, intensity));
                    index += 1;
                    continue;
                }
                let entry_rt = profile.entries[current].0;
                if entry_rt > peak.rt {
                    previous = profile.insert_after(previous, (peak.rt, intensity));
                    index += 1;
                } else if entry_rt < peak.rt {
                    previous = current;
                    current = profile.next[current];
                } else if entry_rt == peak.rt {
                    profile.entries[current].1 = x86_64::add(intensity, profile.entries[current].1);
                    previous = current;
                    current = profile.next[current];
                    index += 1;
                } else {
                    return Err(Error::InvalidValue(NAN_RT_MERGE_WHAT.into()));
                }
            }
        }
        Ok(profile.into_vec())
    }

    /// Total peak count, after checking both profile ceilings.
    fn preflight_profile(&self) -> Result<usize> {
        let mut peaks = 0usize;
        let mut steps = 0usize;
        for (index, trace) in self.traces.iter().enumerate() {
            let count = trace.peaks.len();
            if index > 0 {
                steps = steps.saturating_add(peaks).saturating_add(count);
            }
            peaks = peaks.saturating_add(count);
        }
        preflight_peaks(peaks)?;
        if steps > Self::MAX_PROFILE_STEPS {
            return Err(Error::InvalidValue(format!(
                "an intensity profile over {} traces and {peaks} peaks may take {steps} merge \
                 steps, exceeding the limit of {}",
                self.traces.len(),
                Self::MAX_PROFILE_STEPS
            )));
        }
        Ok(peaks)
    }
}

impl Index<usize> for MassTraces {
    type Output = MassTrace;

    /// The trace at `index`: source `operator[]`. Panics when out of range, like
    /// slice indexing, where the source is undefined; [`MassTraces::get`] is the
    /// checked accessor.
    fn index(&self, index: usize) -> &MassTrace {
        &self.traces[index]
    }
}

impl IndexMut<usize> for MassTraces {
    /// The mutable trace at `index`. Panics when out of range, like slice
    /// indexing; [`MassTraces::get_mut`] is the checked accessor.
    fn index_mut(&mut self, index: usize) -> &mut MassTrace {
        &mut self.traces[index]
    }
}

impl<'a> IntoIterator for &'a MassTraces {
    type Item = &'a MassTrace;
    type IntoIter = std::slice::Iter<'a, MassTrace>;

    fn into_iter(self) -> Self::IntoIter {
        self.traces.iter()
    }
}

impl<'a> IntoIterator for &'a mut MassTraces {
    type Item = &'a mut MassTrace;
    type IntoIter = std::slice::IterMut<'a, MassTrace>;

    fn into_iter(self) -> Self::IntoIter {
        self.traces.iter_mut()
    }
}

/// A theoretical isotope pattern (source
/// `FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern`).
///
/// The source declares no constructor, so a default-initialised pattern holds
/// indeterminate `optional_begin`, `optional_end`, `max` and `trimmed_left`
/// until the algorithm assigns them; they are zero here.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TheoreticalIsotopePattern {
    /// Intensity contribution of each isotope peak, `f64` like the source
    /// `std::vector<double>`.
    pub intensity: Vec<f64>,
    /// Number of optional peaks at the beginning of the pattern.
    pub optional_begin: usize,
    /// Number of optional peaks at the end of the pattern.
    pub optional_end: usize,
    /// The maximum intensity contribution before the pattern was scaled to 1.
    pub max: f64,
    /// Number of isotopes trimmed on the left side, needed to reconstruct the
    /// monoisotopic peak.
    pub trimmed_left: usize,
}

impl TheoreticalIsotopePattern {
    /// The number of isotope peaks, the length of [`Self::intensity`]: source
    /// `size`.
    pub fn len(&self) -> usize {
        self.intensity.len()
    }

    /// Whether the pattern holds no isotope peaks.
    pub fn is_empty(&self) -> bool {
        self.intensity.is_empty()
    }
}

/// The peak matched to one isotope of an [`IsotopePattern`].
///
/// The source stores a signed index: `-1` when no peak was found and `-2` when
/// the peak was removed to improve the isotope fit. It assigns no other
/// negative value, so the three cases map to variants without loss.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum PatternPeak {
    /// No peak was found for the isotope: source `-1`, the initial value.
    #[default]
    NotFound,
    /// The peak was removed to improve the isotope fit: source `-2`.
    Removed,
    /// Index of the matched peak within the spectrum at the same position of
    /// [`IsotopePattern::spectrum`].
    Found(usize),
}

impl PatternPeak {
    /// The peak index when a peak is matched, otherwise `None`.
    pub const fn index(self) -> Option<usize> {
        match self {
            Self::Found(index) => Some(index),
            Self::NotFound | Self::Removed => None,
        }
    }
}

/// An isotope pattern found in the data (source
/// `FeatureFinderAlgorithmPickedHelperStructs::IsotopePattern`).
///
/// The five vectors are indexed by isotope and have equal length after
/// [`Self::new`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IsotopePattern {
    /// Matched peak of each isotope.
    pub peak: Vec<PatternPeak>,
    /// Spectrum index of each isotope; meaningless unless [`Self::peak`] holds
    /// [`PatternPeak::Found`] at the same position, as the source documents.
    pub spectrum: Vec<usize>,
    /// Peak intensity of each isotope; `0` when no peak is matched.
    pub intensity: Vec<f64>,
    /// m/z score of each isotope's peak; `0` when no peak is matched.
    pub mz_score: Vec<f64>,
    /// Theoretical m/z of each isotope peak.
    pub theoretical_mz: Vec<f64>,
    /// The theoretical isotope pattern.
    pub theoretical_pattern: TheoreticalIsotopePattern,
}

impl IsotopePattern {
    /// Most isotopes [`Self::new`] allocates for.
    pub const MAX_SIZE: usize = 100_000;

    /// A pattern of `size` isotopes with no peak matched: source
    /// `IsotopePattern(Size size)`.
    ///
    /// Every [`Self::peak`] is [`PatternPeak::NotFound`], the source's `-1`, and
    /// every spectrum index, intensity, m/z score and theoretical m/z is zero,
    /// as in the source's value-initialised vectors. The
    /// [`Self::theoretical_pattern`] is empty with zero scalars; the source
    /// default-initialises it, leaving its scalars indeterminate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `size` exceeds [`Self::MAX_SIZE`],
    /// checked before allocating.
    pub fn new(size: usize) -> Result<Self> {
        if size > Self::MAX_SIZE {
            return Err(Error::InvalidValue(format!(
                "an isotope pattern of {size} isotopes exceeds the limit of {}",
                Self::MAX_SIZE
            )));
        }
        Ok(Self {
            peak: vec![PatternPeak::NotFound; size],
            spectrum: vec![0; size],
            intensity: vec![0.0; size],
            mz_score: vec![0.0; size],
            theoretical_mz: vec![0.0; size],
            theoretical_pattern: TheoreticalIsotopePattern::default(),
        })
    }
}

/// End marker of [`LinkedProfile`]; never a node index, because node counts are
/// bounded by [`MassTraces::MAX_PEAKS`].
const NIL: usize = usize::MAX;

/// The source's `std::list` profile as an index-linked list.
struct LinkedProfile {
    entries: Vec<(f64, f64)>,
    next: Vec<usize>,
    head: usize,
}

impl LinkedProfile {
    fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
            next: Vec::with_capacity(capacity),
            head: NIL,
        }
    }

    /// Links `entry` after `previous`, or at the head when `previous` is
    /// [`NIL`], and returns the new node.
    fn insert_after(&mut self, previous: usize, entry: (f64, f64)) -> usize {
        let node = self.entries.len();
        self.entries.push(entry);
        if previous == NIL {
            self.next.push(self.head);
            self.head = node;
        } else {
            self.next.push(self.next[previous]);
            self.next[previous] = node;
        }
        node
    }

    fn into_vec(self) -> Vec<(f64, f64)> {
        let mut profile = Vec::with_capacity(self.entries.len());
        let mut node = self.head;
        while node != NIL {
            profile.push(self.entries[node]);
            node = self.next[node];
        }
        profile
    }
}

fn precondition(message: &str) -> Error {
    Error::InvalidValue(message.to_owned())
}

fn preflight_peaks(count: usize) -> Result<()> {
    if count > MassTraces::MAX_PEAKS {
        return Err(Error::InvalidValue(format!(
            "{count} mass-trace peaks exceed the limit of {}",
            MassTraces::MAX_PEAKS
        )));
    }
    Ok(())
}
