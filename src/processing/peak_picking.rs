// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! High-resolution centroiding with histogram noise estimation and natural cubic
//! splines.
//!
//! Ports `src/openms/include/OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h`
//! and `src/openms/source/PROCESSING/CENTROIDING/PeakPickerHiRes.cpp`, with the
//! signal-to-noise estimator of `SignalToNoiseEstimatorMedian.h` in the `noise`
//! submodule. `docs/PEAK_PICKING_SUPPORT.md` lists every source member and its
//! counterpart, the preserved source conventions, the native differences and
//! the evidence.
//!
//! The source class detects ion signals in profile data and reconstructs each
//! peak shape by cubic-spline interpolation; a picked peak's position and
//! intensity are those of the spline maximum. It was tested mainly on
//! high-resolution data (FT-ICR, Orbitrap); with noise reduction and baseline
//! subtraction it may also suit low-resolution data. Picking is
//! one-dimensional: ion-mobility-separated data is treated as one-dimensional
//! data and the picked peaks report the intensity-weighted ion mobility, which
//! is correct for binned data and wrong for fully two-dimensional input.
//!
//! The source notes that peaks must be sorted by ascending position but does
//! not check it. This port checks it by default and returns
//! [`Error::UnsortedData`] otherwise.
//!
//! # Native default and source compatibility
//!
//! The native default refuses inputs on which the source's behaviour is
//! degenerate: unsorted or duplicate positions, negative intensities, a
//! non-positive spline maximum, a ppm width at a non-positive position and a
//! second ion mobility array. [`PickingCompatibility::source`](crate::processing::peak_picking::PickingCompatibility::source)
//! selects the source behaviour for each
//! of these; a TOPP tool reproducing C++ output opts in. The picker is serial,
//! as the source is (it has no OpenMP).

mod noise;
pub use noise::{
    NoiseEstimates, NoiseHistogramRange, NoiseRangeParameters,
    SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME, SignalToNoiseEstimatorMedian,
};
// The natural cubic spline now lives with the rest of the MATH/MISC splines in
// `crate::processing::spline`; this re-export keeps its original path, which the
// peak picker and the retention-time transformations import.
pub use super::spline::CubicSpline2d;

use super::spline::bisection::{DEFAULT_BISECTION_THRESHOLD, MAX_BISECTION_STEPS};
use super::spline::spline_bisection;
use super::{SpectrumFilter, checked_intensity};
use crate::kernel::{
    ChromatogramPeak, DataArray, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, SpectrumType,
    SpectrumTypeQueryLimits,
};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use noise::{flag, flag_entry, float, float_entry, int_entry, integer, string_entry};

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

/// Parameter handler name of the source class, used in parameter diagnostics.
pub const PEAK_PICKER_HI_RES_NAME: &str = "PeakPickerHiRes";

/// The exception text of source `pickExperiment` when a selected spectrum is
/// centroided and the spectrum type is checked.
pub const CENTROIDED_INPUT_MESSAGE: &str =
    "Error: Centroided data provided but profile spectra expected.";

/// A profile sample the picker reads: a position and a float intensity.
trait SignalPoint {
    fn position(&self) -> f64;
    fn intensity(&self) -> f32;
}
impl SignalPoint for Peak1D {
    fn position(&self) -> f64 {
        self.mz
    }
    fn intensity(&self) -> f32 {
        self.intensity
    }
}
impl SignalPoint for ChromatogramPeak {
    fn position(&self) -> f64 {
        self.rt
    }
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

/// Check point count, finiteness, order and the compatibility-dependent
/// refusals, in that order.
fn validate_points<F: Fn(usize) -> f64, G: Fn(usize) -> f64>(
    n: usize,
    position: &F,
    intensity: &G,
    max_points: usize,
    compatibility: &PickingCompatibility,
) -> Result<()> {
    if max_points == 0 || n > max_points {
        return Err(bad("signal arrays differ in length or exceed point limit"));
    }
    for i in 0..n {
        let (x, y) = (position(i), intensity(i));
        if !x.is_finite() || !y.is_finite() {
            return Err(bad(
                "peak picking requires finite coordinates and intensities",
            ));
        }
        if y < 0.0 && !compatibility.allow_negative_intensities {
            return Err(bad(
                "peak picking requires nonnegative intensities; negative intensities need PickingCompatibility::allow_negative_intensities",
            ));
        }
    }
    if !compatibility.allow_unsorted_positions {
        for i in 1..n {
            if position(i - 1) > position(i) {
                return Err(Error::UnsortedData);
            }
        }
    }
    if !compatibility.allow_duplicate_positions {
        for i in 1..n {
            if position(i - 1) == position(i) {
                return Err(bad(
                    "peak picking requires distinct coordinates; duplicates need PickingCompatibility::allow_duplicate_positions",
                ));
            }
        }
    }
    Ok(())
}

/// Which source behaviours the picker and its noise estimator adopt where the
/// native default refuses.
///
/// Every flag is `false` by default, and [`PickingCompatibility::source`] sets
/// them all. The pattern follows `dta::WriteOptions::source`: the guard stays
/// the library default and a tool reproducing the C++ output selects the source
/// explicitly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PickingCompatibility {
    /// Accept equal consecutive positions. The source stores each peak's
    /// support in a `std::map` keyed by position, so a later sample at an equal
    /// position overwrites the stored intensity, while the intensity-weighted
    /// ion mobility still accumulates every added sample; a spacing check with a
    /// zero minimum spacing fails. The port reproduces all three.
    pub allow_duplicate_positions: bool,
    /// Accept decreasing positions. The source documents sorted input as a
    /// precondition but does not check it, and its class test picks three
    /// unsorted tandem spectra: the support map sorts each peak's samples by
    /// position, the spacing tests work with whatever signed differences
    /// arise, and a reversed bisection bracket runs one step. The port
    /// reproduces that arithmetic. Duplicate detection under
    /// `allow_duplicate_positions = false` then only sees consecutive samples.
    pub allow_unsorted_positions: bool,
    /// Accept negative intensities, for example baseline-corrected data. The
    /// source picks them like any other value, and its noise estimator puts
    /// them in the first bin or, when the automatic histogram range comes out
    /// negative, reports a signal-to-noise ratio of zero everywhere.
    pub allow_negative_intensities: bool,
    /// Report a centroid whose spline maximum is zero or negative, as the source
    /// does. With FWHM reporting the source then never leaves its half-height
    /// bisection; the port returns an error instead.
    pub allow_nonpositive_maximum: bool,
    /// Report a ppm FWHM at a zero or negative centroid position: the source
    /// divides by the position regardless. A non-finite result is still an
    /// error.
    pub allow_nonpositive_fwhm_position: bool,
    /// Take the first ion mobility float array and give the output array only
    /// its name, as source `pick` does. The native default refuses a spectrum
    /// with more than one matching array and copies the input array's
    /// description (metadata and processing handles) onto the output array.
    pub source_mobility_arrays: bool,
}
impl PickingCompatibility {
    /// Every source behaviour.
    pub const fn source() -> Self {
        Self {
            allow_duplicate_positions: true,
            allow_unsorted_positions: true,
            allow_negative_intensities: true,
            allow_nonpositive_maximum: true,
            allow_nonpositive_fwhm_position: true,
            source_mobility_arrays: true,
        }
    }
}

/// Positions of the first and last profile samples an extended peak reached.
///
/// Source `PeakPickerHiRes::PeakBoundary`, whose members are named `mz_min` and
/// `mz_max` for spectra and carry retention times for chromatograms. The
/// boundary includes a final sample the extension visited and rejected as
/// missing, so it is not necessarily a spline knot or a half-height crossing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PeakBoundary {
    /// Position of the leftmost visited sample (source `mz_min`), in Th for
    /// spectra and seconds for chromatograms.
    pub min: f64,
    /// Position of the rightmost visited sample (source `mz_max`).
    pub max: f64,
}

/// Unit of the reported full width at half maximum.
///
/// Source parameter `report_FWHM_unit`: `absolute` is [`FwhmUnit::Absolute`]
/// and every other value, `relative` in practice, is [`FwhmUnit::Ppm`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FwhmUnit {
    /// Width in the unit of the input position, written to a float array named
    /// `FWHM`.
    Absolute,
    /// Width divided by the centroid position times one million, written to a
    /// float array named `FWHM_ppm`. Only sensible for spectra.
    Ppm,
}

/// A picked spectrum with its peak boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct PickedSpectrum {
    /// The centroided spectrum.
    pub spectrum: MSSpectrum,
    /// One boundary per centroid, in centroid order.
    pub boundaries: Vec<PeakBoundary>,
    /// Names of profile annotation arrays that were dropped because they have
    /// no defined centroid aggregation rule (native report; the source drops
    /// them silently).
    pub omitted_arrays: Vec<String>,
}

/// A picked chromatogram with its peak boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct PickedChromatogram {
    /// The centroided chromatogram.
    pub chromatogram: MSChromatogram,
    /// One boundary per centroid, in centroid order, in seconds.
    pub boundaries: Vec<PeakBoundary>,
    /// Names of dropped profile annotation arrays (native report).
    pub omitted_arrays: Vec<String>,
}

/// The per-record reports of an experiment pick, without the picked experiment.
///
/// Returned by [`PeakPickerHiRes::pick_experiment_in_place`], which centroids the
/// caller's experiment instead of building a second one, so it has no experiment
/// of its own to hand back. The four vectors are those of
/// [`PickedExperiment`] and carry the same per-record meaning.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PickedExperimentReport {
    /// One entry per spectrum; `None` means the spectrum was left as it was.
    pub spectrum_boundaries: Vec<Option<Vec<PeakBoundary>>>,
    /// One entry per chromatogram.
    pub chromatogram_boundaries: Vec<Vec<PeakBoundary>>,
    /// Dropped annotation array names per spectrum.
    pub omitted_spectrum_arrays: Vec<Vec<String>>,
    /// Dropped annotation array names per chromatogram.
    pub omitted_chromatogram_arrays: Vec<Vec<String>>,
}

/// A picked experiment with per-record boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct PickedExperiment {
    /// The output experiment: experimental settings copied, selected spectra
    /// picked, other spectra copied, every chromatogram picked.
    pub experiment: MSExperiment,
    /// One entry per input spectrum; `None` means copied without picking. The
    /// source appends boundaries for picked spectra only.
    pub spectrum_boundaries: Vec<Option<Vec<PeakBoundary>>>,
    /// One entry per chromatogram.
    pub chromatogram_boundaries: Vec<Vec<PeakBoundary>>,
    /// Dropped annotation array names per input spectrum.
    pub omitted_spectrum_arrays: Vec<Vec<String>>,
    /// Dropped annotation array names per chromatogram.
    pub omitted_chromatogram_arrays: Vec<Vec<String>>,
}

/// The high-resolution peak picker, source `PeakPickerHiRes`.
///
/// The typed fields are the source members that `updateMembers_` derives from
/// the parameters, plus native options. [`PeakPickerHiRes::defaults`],
/// [`PeakPickerHiRes::from_param`] and [`PeakPickerHiRes::to_param`] carry the
/// source `DefaultParamHandler` contract.
///
/// # Algorithm
///
/// A peak core is a sample `i` with `2 <= i < n - 2` that is a strict intensity
/// maximum, whose immediate neighbours are nonzero (magnitude at least
/// `f64::EPSILON`), and whose three samples reach the signal-to-noise threshold.
/// With spacing checks, `min_spacing` is the smaller of the two apex spacings,
/// and a neighbour counts as present when its spacing is below
/// `spacing_difference * min_spacing`; both neighbours must be present, or one
/// with [`allow_missing_flank`](Self::allow_missing_flank). A core flanked by
/// two more intense satellites (an oscillation) is skipped together with the
/// next sample. The core then extends to each side along non-increasing
/// intensities, stopping after a zero intensity, at a spacing beyond
/// `spacing_difference_gap * min_spacing`, or once more than
/// [`missing`](Self::missing) samples failed the signal-to-noise or spacing test.
/// At least three support samples are needed. The centroid is the maximum of
/// the natural cubic spline through the support, found by
/// [`spline_bisection`] between the two neighbours (or the core where a
/// neighbour is missing) with a bracket width of `1e-6`. Picking resumes after
/// the last sample the right extension visited.
#[derive(Clone, Debug, PartialEq)]
pub struct PeakPickerHiRes {
    /// Minimal signal-to-noise ratio of a picked peak's core and support,
    /// `signal_to_noise`; `0.0` (the default) disables noise estimation.
    pub signal_to_noise: f64,
    /// The estimator behind `signal_to_noise`, the `SignalToNoise:` subsection.
    pub noise_estimator: SignalToNoiseEstimatorMedian,
    /// `spacing_difference`, default `1.5`: the maximum spacing during
    /// extension in multiples of `min_spacing` before a sample counts as
    /// missing. `0.0` disables the constraint, as the source's mapping to
    /// infinity does. Not applicable to chromatograms.
    pub spacing_difference: f64,
    /// `spacing_difference_gap`, default `4.0`: extension stops when the spacing
    /// to the next sample exceeds this multiple of `min_spacing`. `0.0` disables
    /// the constraint; disabling both constraints disables spacing checks.
    pub spacing_difference_gap: f64,
    /// `missing`, default `1`: missing samples allowed per extension side.
    pub missing: usize,
    /// `allow_missing_flank`, default `false`: accept a core with only one
    /// present neighbour (for TimsTOF data missing a leading or trailing edge).
    pub allow_missing_flank: bool,
    /// `report_FWHM` and `report_FWHM_unit`: `Some` adds a float array with the
    /// FWHM of every centroid.
    pub report_fwhm: Option<FwhmUnit>,
    /// The `report_FWHM_unit` value kept while [`report_fwhm`](Self::report_fwhm)
    /// is `None`, so that [`PeakPickerHiRes::to_param`] reproduces it. Ignored
    /// while `report_fwhm` is `Some`; [`PeakPickerHiRes::from_param`] then
    /// leaves it at its default, [`FwhmUnit::Ppm`] (`relative`).
    pub inactive_fwhm_unit: FwhmUnit,
    /// `ms_levels`: empty (the default) selects automatic mode, which picks
    /// every spectrum not already centroided; otherwise only these MS levels
    /// are picked and others are copied.
    pub ms_levels: Vec<u32>,
    /// The `check_spectrum_type` argument of source `pickExperiment` (not a
    /// parameter): in manual mode, refuse a selected centroided spectrum.
    pub check_spectrum_type: bool,
    /// Native: an exact float array name whose values are intensity-weighted
    /// per centroid. `None` (the default) applies the source rule, the first
    /// array for which [`MSSpectrum::contains_im_data`] holds.
    pub ion_mobility_array: Option<String>,
    /// Native: source behaviours to adopt; see [`PickingCompatibility`].
    pub compatibility: PickingCompatibility,
    /// Native ceiling on the points of one record.
    pub max_points: usize,
    /// Native bound on apex candidates and extension samples for one record.
    pub max_work: usize,
    /// Native per-record allowance of the acquisition-metadata copy ledger, in
    /// both work units and bytes.
    ///
    /// [`PeakPickerHiRes::pick_experiment`] and
    /// [`PeakPickerHiRes::pick_experiment_in_place`] charge every record's
    /// acquisition metadata to one ledger before and while copying it. That
    /// ledger's fixed part alone is a ceiling on the *number* of records rather
    /// than on any one record: a Q Exactive run whose spectra carry the usual
    /// scan window, acquisition and source-file metadata spends about 8 KiB of
    /// it per spectrum, so the fixed part alone stops at roughly 34 000 spectra,
    /// well inside the size of an ordinary LC-MS run and with no counterpart in
    /// source `pickExperiment`. This allowance is added to the fixed part once
    /// per input record, which makes the ledger track the input instead of
    /// capping it, while a record whose metadata dwarfs the whole input is still
    /// refused.
    ///
    /// The default, 64 KiB, is about eight times the per-spectrum cost of an
    /// ordinary vendor-converted run. Zero pins the ledger at its fixed part,
    /// which is the behaviour before this field existed.
    ///
    /// Like [`max_points`](Self::max_points) and [`max_work`](Self::max_work)
    /// this is a Rust-API-only field: it is not in
    /// [`PeakPickerHiRes::defaults`], so [`PeakPickerHiRes::to_param`] does not
    /// emit it and [`PeakPickerHiRes::from_param`] cannot set it, and a caller
    /// driving the picker from a TOPP `.ini` therefore always gets the default.
    /// A converter whose per-spectrum metadata is richer than the run the
    /// default is calibrated on has to raise it through the Rust API.
    pub max_metadata_per_record: usize,
}
impl Default for PeakPickerHiRes {
    fn default() -> Self {
        Self {
            signal_to_noise: 0.0,
            noise_estimator: Default::default(),
            spacing_difference: 1.5,
            spacing_difference_gap: 4.0,
            missing: 1,
            allow_missing_flank: false,
            report_fwhm: None,
            inactive_fwhm_unit: FwhmUnit::Ppm,
            ms_levels: Vec::new(),
            check_spectrum_type: true,
            ion_mobility_array: None,
            compatibility: PickingCompatibility::default(),
            max_points: 1_000_000,
            max_work: 10_000_000,
            max_metadata_per_record: 64 * 1024,
        }
    }
}

struct PickedSignal {
    positions: Vec<f64>,
    intensities: Vec<f32>,
    boundaries: Vec<PeakBoundary>,
    fwhm: Vec<f32>,
    mobility: Vec<f32>,
}

/// The source's per-peak `std::map<double, double>` of support samples, kept in
/// two reusable vectors so that picking a peak allocates nothing once they have
/// grown.
///
/// Invariant: `left` holds keys below `right[0]` in descending order, `right`
/// holds the core and larger keys in ascending order, so the map order is
/// `left` reversed followed by `right`. For position-ordered input the core and
/// the right extension only append at or above the largest key and the left
/// extension only at or below the smallest, which are the O(1) paths. A key
/// equal to a stored key overwrites that value and keeps the stored key, as
/// `map[key] = value` does (keys compare with `<`, so `-0.0` equals `0.0`).
/// Only unsorted input reaches the binary-search insertion.
#[derive(Default)]
struct Support {
    /// Keys below `right[0]`, descending.
    left: Vec<(f64, f64)>,
    /// The core and the keys above it, ascending.
    right: Vec<(f64, f64)>,
    xs: Vec<f64>,
    ys: Vec<f64>,
}
impl Support {
    fn reset(&mut self, key: f64, value: f64) {
        self.left.clear();
        self.right.clear();
        self.right.push((key, value));
    }
    fn first_mut(&mut self) -> Option<&mut (f64, f64)> {
        if self.left.is_empty() {
            self.right.first_mut()
        } else {
            self.left.last_mut()
        }
    }
    fn first(&self) -> (f64, f64) {
        match (self.left.last(), self.right.first()) {
            (Some(&pair), _) | (None, Some(&pair)) => pair,
            (None, None) => (f64::NAN, f64::NAN),
        }
    }
    fn last(&self) -> (f64, f64) {
        self.right.last().copied().unwrap_or((f64::NAN, f64::NAN))
    }
    fn insert(&mut self, key: f64, value: f64) {
        let (first_key, _) = self.first();
        let (last_key, _) = self.last();
        if key == first_key {
            if let Some(pair) = self.first_mut() {
                pair.1 = value;
            }
        } else if key < first_key {
            self.left.push((key, value));
        } else if key == last_key {
            if let Some(pair) = self.right.last_mut() {
                pair.1 = value;
            }
        } else if key > last_key {
            self.right.push((key, value));
        } else {
            // Strictly between the extremes: unsorted input only. Keys are
            // finite, so partial_cmp always orders them.
            let order = |a: f64, b: f64| a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal);
            let core = self.right.first().map_or(f64::NAN, |pair| pair.0);
            if key >= core {
                match self.right.binary_search_by(|pair| order(pair.0, key)) {
                    Ok(i) => self.right[i].1 = value,
                    Err(i) => self.right.insert(i, (key, value)),
                }
            } else {
                match self.left.binary_search_by(|pair| order(key, pair.0)) {
                    Ok(i) => self.left[i].1 = value,
                    Err(i) => self.left.insert(i, (key, value)),
                }
            }
        }
    }
    fn len(&self) -> usize {
        self.left.len() + self.right.len()
    }
    /// Fill `xs`/`ys` in key order.
    fn flatten(&mut self) {
        self.xs.clear();
        self.ys.clear();
        for &(key, value) in self.left.iter().rev().chain(self.right.iter()) {
            self.xs.push(key);
            self.ys.push(value);
        }
    }
}

/// One side of the source's half-height bisection, which the source writes
/// with a different midpoint expression per side.
#[derive(Clone, Copy)]
enum HalfHeightSide {
    /// `mz_left / 2 + mz_center / 2`.
    Left,
    /// `(mz_right + mz_center) / 2`.
    Right,
}

fn half_height(
    spline: &CubicSpline2d,
    mut edge: f64,
    mut center: f64,
    height: f64,
    threshold: f64,
    side: HalfHeightSide,
) -> Result<f64> {
    for _ in 0..MAX_BISECTION_STEPS {
        let mid = match side {
            HalfHeightSide::Left => edge / 2.0 + center / 2.0,
            HalfHeightSide::Right => (edge + center) / 2.0,
        };
        let value = spline.eval(mid)?;
        // do { ... } while (fabs(int_mid - fwhm_int) > threshold); the spline
        // evaluation is finite, so the negation is a plain `<=`.
        if (value - height).abs() <= threshold {
            return Ok(mid);
        }
        // Once the midpoint equals a bracket end it never moves again, and the
        // source's do/while never terminates.
        if mid == edge || mid == center {
            return Err(bad(
                "FWHM bisection cannot reach its intensity tolerance (the source does not terminate here)",
            ));
        }
        if value < height {
            edge = mid;
        } else {
            center = mid;
        }
    }
    Err(bad("FWHM bisection exceeded its iteration ceiling"))
}

impl PeakPickerHiRes {
    fn validate(&self) -> Result<()> {
        let nonnegative = |v: f64| !v.is_nan() && v >= 0.0;
        if !self.signal_to_noise.is_finite()
            || self.signal_to_noise < 0.0
            || !nonnegative(self.spacing_difference)
            || !nonnegative(self.spacing_difference_gap)
            || self.max_points == 0
            || self.max_work == 0
            || self.ms_levels.contains(&0)
        {
            return Err(bad("invalid peak picker parameters or resource limits"));
        }
        Ok(())
    }

    /// Source `pick_`, over the samples of one record.
    fn pick_signal<P: SignalPoint>(
        &self,
        points: &[P],
        mut check_spacings: bool,
        mobility: Option<&[f32]>,
    ) -> Result<PickedSignal> {
        self.validate()?;
        let n = points.len();
        let x = |j: usize| points[j].position();
        let y32 = |j: usize| points[j].intensity();
        let y = |j: usize| f64::from(points[j].intensity());
        validate_points(n, &x, &y, self.max_points, &self.compatibility)?;
        if let Some(values) = mobility {
            if values.len() != n || values.iter().any(|v| !v.is_finite()) {
                return Err(bad(
                    "ion mobility must be finite and aligned with input peaks",
                ));
            }
        }
        let mut out = PickedSignal {
            positions: Vec::new(),
            intensities: Vec::new(),
            boundaries: Vec::new(),
            fwhm: Vec::new(),
            mobility: Vec::new(),
        };
        if n < 5 {
            return Ok(out);
        }
        // updateMembers_: zero disables a constraint by mapping it to infinity.
        let spacing = if self.spacing_difference == 0.0 {
            f64::INFINITY
        } else {
            self.spacing_difference
        };
        let gap = if self.spacing_difference_gap == 0.0 {
            f64::INFINITY
        } else {
            self.spacing_difference_gap
        };
        if spacing == f64::INFINITY && gap == f64::INFINITY {
            check_spacings = false;
        }
        let estimates = if self.signal_to_noise > 0.0 {
            Some(
                self.noise_estimator
                    .estimate_signal(points, &self.compatibility)?,
            )
        } else {
            None
        };
        // With estimation disabled the source compares 0.0 >= 0.0.
        let passes_sn = |j: usize| {
            estimates
                .as_ref()
                .is_none_or(|s| s.signal_to_noise[j] >= self.signal_to_noise)
        };
        let weight = |j: usize, values: &[f32]| f64::from(values[j] * y32(j));
        let mut used = 0usize;
        let mut charge = || -> Result<()> {
            if used == self.max_work {
                Err(bad("peak picking exceeds configured work limit"))
            } else {
                used += 1;
                Ok(())
            }
        };
        let mut support = Support::default();
        let mut i = 2;
        while i < n - 2 {
            charge()?;
            let (central_x, central_y) = (x(i), y(i));
            let (left_x, left_y) = (x(i - 1), y(i - 1));
            let (right_x, right_y) = (x(i + 1), y(i + 1));
            // Do not interpolate when a neighbour is a zero data point.
            if left_y.abs() < f64::EPSILON || right_y.abs() < f64::EPSILON {
                i += 1;
                continue;
            }
            let (mut has_left, mut has_right, mut min_spacing) = (false, false, 0.0);
            if check_spacings {
                let left_to_central = central_x - left_x;
                let central_to_right = right_x - central_x;
                min_spacing = if left_to_central < central_to_right {
                    left_to_central
                } else {
                    central_to_right
                };
                has_left = left_to_central < spacing * min_spacing;
                has_right = central_to_right < spacing * min_spacing;
            }
            let spacing_ok = !check_spacings
                || if self.allow_missing_flank {
                    has_left || has_right
                } else {
                    has_left && has_right
                };
            if !(central_y > left_y
                && central_y > right_y
                && passes_sn(i)
                && passes_sn(i - 1)
                && passes_sn(i + 1)
                && spacing_ok)
            {
                i += 1;
                continue;
            }
            // A core surrounded by more intense satellites is an oscillation.
            let has_left2 = left_x - x(i - 2) < spacing * min_spacing;
            let has_right2 = x(i + 2) - right_x < spacing * min_spacing;
            let spacing_ok2 = !check_spacings
                || if self.allow_missing_flank {
                    has_left2 || has_right2
                } else {
                    has_left2 && has_right2
                };
            if left_y < y(i - 2)
                && right_y < y(i + 2)
                && passes_sn(i - 2)
                && passes_sn(i + 2)
                && spacing_ok2
            {
                i += 2;
                continue;
            }
            let include_left = !check_spacings || has_left;
            let include_right = !check_spacings || has_right;
            support.reset(central_x, central_y);
            let mut weighted_im = 0.0;
            if let Some(values) = mobility {
                weighted_im += weight(i, values);
            }
            if include_left {
                support.insert(left_x, left_y);
            }
            if include_right {
                support.insert(right_x, right_y);
            }
            if let Some(values) = mobility {
                if include_left {
                    weighted_im += weight(i - 1, values);
                }
                if include_right {
                    weighted_im += weight(i + 1, values);
                }
            }

            // Extend to the left.
            let mut k = 2;
            let mut previous_zero = false;
            let mut missing = 0usize;
            let mut left_boundary = i - 1;
            while k <= i && !previous_zero && missing <= self.missing {
                let j = i - k;
                let (first_x, first_y) = support.first();
                if !(y(j) <= first_y && (!check_spacings || first_x - x(j) < gap * min_spacing)) {
                    break;
                }
                charge()?;
                let good =
                    passes_sn(j) && (!check_spacings || first_x - x(j) < spacing * min_spacing);
                if !good {
                    missing = missing
                        .checked_add(1)
                        .ok_or_else(|| bad("missing-point count overflow"))?;
                }
                if good || missing <= self.missing {
                    support.insert(x(j), y(j));
                    if let Some(values) = mobility {
                        weighted_im += weight(j, values);
                    }
                }
                previous_zero = y(j) == 0.0;
                left_boundary = j;
                k += 1;
            }

            // Extend to the right.
            k = 2;
            previous_zero = false;
            missing = 0;
            let mut right_boundary = i + 1;
            while i + k < n && !previous_zero && missing <= self.missing {
                let j = i + k;
                let (last_x, last_y) = support.last();
                if !(y(j) <= last_y && (!check_spacings || x(j) - last_x < gap * min_spacing)) {
                    break;
                }
                charge()?;
                let good =
                    passes_sn(j) && (!check_spacings || x(j) - last_x < spacing * min_spacing);
                if !good {
                    missing = missing
                        .checked_add(1)
                        .ok_or_else(|| bad("missing-point count overflow"))?;
                }
                if good || missing <= self.missing {
                    support.insert(x(j), y(j));
                    if let Some(values) = mobility {
                        weighted_im += weight(j, values);
                    }
                }
                previous_zero = y(j) == 0.0;
                right_boundary = j;
                k += 1;
            }

            if support.len() < 3 {
                i += 1;
                continue;
            }
            support.flatten();
            let spline = CubicSpline2d::with_max_points(&support.xs, &support.ys, self.max_points)?;
            // A missing neighbour brackets the maximum at the core instead.
            let (bracket_left, bracket_right) = if check_spacings {
                (
                    if has_left { left_x } else { central_x },
                    if has_right { right_x } else { central_x },
                )
            } else {
                (left_x, right_x)
            };
            let (position, intensity) = spline_bisection(
                &spline,
                bracket_left,
                bracket_right,
                DEFAULT_BISECTION_THRESHOLD,
            )?;
            if intensity <= 0.0 && !self.compatibility.allow_nonpositive_maximum {
                return Err(bad(
                    "peak spline maximum is not positive; this needs PickingCompatibility::allow_nonpositive_maximum",
                ));
            }
            if let Some(unit) = self.report_fwhm {
                let height = intensity / 2.0;
                let threshold = 0.01 * height;
                let lowest = support.xs[0];
                let highest = support.xs[support.xs.len() - 1];
                // The support ending above half height gives its end as the
                // crossing, probably underestimating the width.
                let left = if spline.eval(lowest)? > height {
                    lowest
                } else {
                    half_height(
                        &spline,
                        lowest,
                        position,
                        height,
                        threshold,
                        HalfHeightSide::Left,
                    )?
                };
                let right = if spline.eval(highest)? > height {
                    highest
                } else {
                    half_height(
                        &spline,
                        highest,
                        position,
                        height,
                        threshold,
                        HalfHeightSide::Right,
                    )?
                };
                let width = right - left;
                let width = match unit {
                    FwhmUnit::Absolute => width,
                    FwhmUnit::Ppm => {
                        if position <= 0.0 && !self.compatibility.allow_nonpositive_fwhm_position {
                            return Err(bad(
                                "ppm FWHM needs a positive centroid position; this needs PickingCompatibility::allow_nonpositive_fwhm_position",
                            ));
                        }
                        width / position * 1e6
                    }
                };
                out.fwhm.push(checked_intensity(width)?);
            }
            if mobility.is_some() {
                let mut total = 0.0;
                for &value in &support.ys {
                    total += value;
                }
                out.mobility.push(checked_intensity(weighted_im / total)?);
            }
            out.positions.push(position);
            out.intensities.push(checked_intensity(intensity)?);
            out.boundaries.push(PeakBoundary {
                min: x(left_boundary),
                max: x(right_boundary),
            });
            // Jump over the samples the right extension visited.
            i += k;
        }
        Ok(out)
    }

    /// The float data array whose values are intensity-weighted, if any.
    fn mobility_array(&self, input: &MSSpectrum) -> Result<Option<usize>> {
        let first_only = self.compatibility.source_mobility_arrays;
        if let Some(name) = &self.ion_mobility_array {
            let mut found = None;
            for (index, array) in input.float_data_arrays.iter().enumerate() {
                if array.name == *name {
                    if found.is_none() {
                        found = Some(index);
                        if first_only {
                            break;
                        }
                    } else {
                        return Err(bad("multiple matching ion mobility arrays"));
                    }
                }
            }
            return found
                .map(Some)
                .ok_or_else(|| bad("requested ion mobility array is missing"));
        }
        if !input.contains_im_data() {
            return Ok(None);
        }
        let (index, _) = input.im_data()?;
        if !first_only {
            // `contains_im_data` inspects names only, so a name-only probe
            // applies the same rule to each later array.
            for array in &input.float_data_arrays[index + 1..] {
                let probe = MSSpectrum {
                    float_data_arrays: vec![DataArray::new(array.name.clone(), Vec::new())],
                    ..MSSpectrum::default()
                };
                if probe.contains_im_data() {
                    return Err(bad(
                        "multiple matching ion mobility arrays; the first is used with PickingCompatibility::source_mobility_arrays",
                    ));
                }
            }
        }
        Ok(Some(index))
    }

    /// Pick a profile spectrum with source spacing checks enabled.
    ///
    /// Source `pick(const MSSpectrum&, MSSpectrum&, std::vector<PeakBoundary>&)`,
    /// which also serves the overload without boundaries.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerHiRes::pick_spectrum_with_spacing`].
    pub fn pick_spectrum(&self, input: &MSSpectrum) -> Result<PickedSpectrum> {
        self.pick_spectrum_with_spacing(input, true)
    }

    /// Pick a profile spectrum, choosing whether spacing constraints apply.
    ///
    /// Source `pick(input, output, boundaries, check_spacings)`. The output
    /// keeps every metadata field of the input (source `copySpectrumMeta`),
    /// is marked centroided, and holds only the picked float arrays: the
    /// intensity-weighted ion mobility first, when the input has an ion mobility
    /// array, then `FWHM` or `FWHM_ppm` when reported. Both arrays exist, empty,
    /// even for a spectrum with fewer than five samples, as in the source.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] for invalid options or limits, an invalid
    ///   spectrum, non-finite values, a refusal that [`PickingCompatibility`]
    ///   can lift, a missing, duplicated or misaligned ion mobility array, an
    ///   exhausted work budget, a spline or bisection failure, a non-finite
    ///   narrowed output value, or an FWHM search that cannot converge (the
    ///   source loops forever there).
    /// * [`Error::UnsortedData`] for decreasing positions, unless
    ///   [`PickingCompatibility::allow_unsorted_positions`] is set.
    /// * [`Error::Unsupported`] when noise estimation runs with
    ///   [`NoiseHistogramRange::Percentile`].
    pub fn pick_spectrum_with_spacing(
        &self,
        input: &MSSpectrum,
        check_spacings: bool,
    ) -> Result<PickedSpectrum> {
        self.pick_spectrum_with_acquisition(
            input,
            check_spacings,
            &mut super::AcquisitionCopies::default(),
        )
    }

    pub(super) fn pick_spectrum_with_acquisition(
        &self,
        input: &MSSpectrum,
        check_spacings: bool,
        copies: &mut super::AcquisitionCopies,
    ) -> Result<PickedSpectrum> {
        self.validate()?;
        if input.len() > self.max_points {
            return Err(bad("spectrum exceeds peak picker point limit"));
        }
        input.validate()?;
        let mobility = self.mobility_array(input)?;
        let picked = self.pick_signal(
            &input.peaks,
            check_spacings,
            mobility.map(|j| input.float_data_arrays[j].data.as_slice()),
        )?;
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .enumerate()
            .filter(|(j, _)| Some(*j) != mobility)
            .map(|(_, a)| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        copies.spectrum(input)?;
        let mut output = input.clone();
        output.peaks = picked
            .positions
            .iter()
            .zip(&picked.intensities)
            .map(|(&x, &y)| Peak1D::new(x, y))
            .collect();
        output.spectrum_type = SpectrumType::Centroid;
        output.float_data_arrays.clear();
        output.integer_data_arrays.clear();
        output.string_data_arrays.clear();
        if let Some(j) = mobility {
            let mut array =
                DataArray::new(input.float_data_arrays[j].name.clone(), picked.mobility);
            if !self.compatibility.source_mobility_arrays {
                input.float_data_arrays[j].copy_description_to(&mut array);
            }
            output.float_data_arrays.push(array);
        }
        if let Some(unit) = self.report_fwhm {
            output
                .float_data_arrays
                .push(DataArray::new(fwhm_name(unit), picked.fwhm));
        }
        Ok(PickedSpectrum {
            spectrum: output,
            boundaries: picked.boundaries,
            omitted_arrays,
        })
    }

    /// Pick a chromatogram with spacing checks disabled, as the source overload
    /// defaults.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerHiRes::pick_chromatogram_with_spacing`].
    pub fn pick_chromatogram(&self, input: &MSChromatogram) -> Result<PickedChromatogram> {
        self.pick_chromatogram_with_spacing(input, false)
    }

    /// Pick a chromatogram, choosing whether spacing constraints apply.
    ///
    /// Source `pick(const MSChromatogram&, MSChromatogram&, boundaries,
    /// check_spacings)`. The output keeps the chromatogram settings, metadata
    /// and name, and holds only the `FWHM`/`FWHM_ppm` array when reported.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerHiRes::pick_spectrum_with_spacing`], without the ion
    /// mobility cases.
    pub fn pick_chromatogram_with_spacing(
        &self,
        input: &MSChromatogram,
        check_spacings: bool,
    ) -> Result<PickedChromatogram> {
        self.pick_chromatogram_with_acquisition(
            input,
            check_spacings,
            &mut super::AcquisitionCopies::default(),
        )
    }

    pub(super) fn pick_chromatogram_with_acquisition(
        &self,
        input: &MSChromatogram,
        check_spacings: bool,
        copies: &mut super::AcquisitionCopies,
    ) -> Result<PickedChromatogram> {
        self.validate()?;
        if input.len() > self.max_points {
            return Err(bad("chromatogram exceeds peak picker point limit"));
        }
        input.validate()?;
        let picked = self.pick_signal(&input.peaks, check_spacings, None)?;
        let omitted_arrays = input
            .float_data_arrays
            .iter()
            .map(|a| a.name.clone())
            .chain(input.integer_data_arrays.iter().map(|a| a.name.clone()))
            .chain(input.string_data_arrays.iter().map(|a| a.name.clone()))
            .collect();
        copies.chromatogram(input)?;
        let mut output = input.clone();
        output.peaks = picked
            .positions
            .iter()
            .zip(&picked.intensities)
            .map(|(&x, &y)| ChromatogramPeak::new(x, y))
            .collect();
        output.float_data_arrays.clear();
        output.integer_data_arrays.clear();
        output.string_data_arrays.clear();
        if let Some(unit) = self.report_fwhm {
            output
                .float_data_arrays
                .push(DataArray::new(fwhm_name(unit), picked.fwhm));
        }
        Ok(PickedChromatogram {
            chromatogram: output,
            boundaries: picked.boundaries,
            omitted_arrays,
        })
    }

    /// Pick all chromatograms and the selected spectra without mutating the
    /// input.
    ///
    /// Source `pickExperiment(input, output, boundaries_spec, boundaries_chrom,
    /// check_spectrum_type)`, with [`check_spectrum_type`](Self::check_spectrum_type)
    /// as the last argument. Spectra are visited in order:
    ///
    /// * **Automatic mode** (empty [`ms_levels`](Self::ms_levels)): a spectrum
    ///   whose [`MSSpectrum::get_type`] with data inspection is centroid is
    ///   copied; every other spectrum is picked. That query takes the stored
    ///   type first, then a `PeakPicking` processing record, then the
    ///   `PeakTypeEstimator` heuristic.
    /// * **Manual mode**: a spectrum of an unlisted MS level is copied. A listed
    ///   spectrum whose type query gives centroid is refused when
    ///   `check_spectrum_type` holds, and picked otherwise.
    ///
    /// Every chromatogram is then picked without spacing checks. The source logs
    /// the picked and total spectra per MS level; this port does not log, and
    /// the boundary entries carry the same information.
    ///
    /// The output is built record by record, as the source's is, so the call
    /// holds the input and the centroids it has produced and never a second copy
    /// of the profile data. What it does still own is a copy of every record it
    /// does not pick, which is what returning an owned experiment from a
    /// borrowed one means; [`PeakPickerHiRes::pick_experiment_in_place`] is the
    /// streaming entry point that avoids even that, at the cost of atomicity.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] with [`CENTROIDED_INPUT_MESSAGE`] for a refused
    ///   centroided spectrum (the source's `Exception::IllegalArgument`).
    /// * Any error of the spectrum or chromatogram picking, or of the type
    ///   query, including its resource limits. An error leaves no partial
    ///   result.
    pub fn pick_experiment(&self, input: &MSExperiment) -> Result<PickedExperiment> {
        let limits = self.type_query_limits();
        let mut copies = self.start_experiment(input)?;
        let mut result = PickedExperiment {
            // Source `pickExperiment` copies the experimental settings, resizes
            // the output to the input and then fills record by record; it never
            // holds a second copy of the profile data. Building the output the
            // same way, instead of cloning the input and overwriting each picked
            // record, is what keeps the peak memory of a multi-gigabyte run at
            // the size of the input plus its centroids.
            experiment: MSExperiment {
                spectra: Vec::with_capacity(input.spectra.len()),
                chromatograms: Vec::with_capacity(input.chromatograms.len()),
                settings: input.settings.clone(),
                sql_run_id: input.sql_run_id,
            },
            spectrum_boundaries: Vec::with_capacity(input.spectra.len()),
            chromatogram_boundaries: Vec::with_capacity(input.chromatograms.len()),
            omitted_spectrum_arrays: Vec::with_capacity(input.spectra.len()),
            omitted_chromatogram_arrays: Vec::with_capacity(input.chromatograms.len()),
        };
        for spectrum in &input.spectra {
            if !self.selects(spectrum, limits)? {
                // Source `output[scan_idx] = input[scan_idx]` for a record that
                // is not picked.
                result.experiment.spectra.push(spectrum.clone());
                result.spectrum_boundaries.push(None);
                result.omitted_spectrum_arrays.push(Vec::new());
                continue;
            }
            let picked = self.pick_spectrum_with_acquisition(spectrum, true, &mut copies)?;
            result.experiment.spectra.push(picked.spectrum);
            result.spectrum_boundaries.push(Some(picked.boundaries));
            result.omitted_spectrum_arrays.push(picked.omitted_arrays);
        }
        for chromatogram in &input.chromatograms {
            let picked =
                self.pick_chromatogram_with_acquisition(chromatogram, false, &mut copies)?;
            result.experiment.chromatograms.push(picked.chromatogram);
            result.chromatogram_boundaries.push(picked.boundaries);
            result
                .omitted_chromatogram_arrays
                .push(picked.omitted_arrays);
        }
        Ok(result)
    }

    /// Centroid an experiment in place, replacing every record as it is picked.
    ///
    /// The streaming form of [`PeakPickerHiRes::pick_experiment`]: it picks the
    /// same records, in the same order, by the same rules, and produces
    /// bit-identical centroids, but it writes each picked record back over its
    /// own profile record instead of collecting a second experiment. The profile
    /// samples of a spectrum are therefore released as soon as its centroids
    /// exist, which is what an mzML-to-mzML tool wants: the peak memory is that
    /// of the loaded run, falling towards the size of its centroids, where
    /// [`PeakPickerHiRes::pick_experiment`] additionally holds the output.
    ///
    /// This entry point exists because the borrowed-input signature of
    /// [`PeakPickerHiRes::pick_experiment`] cannot avoid owning a second copy of
    /// every record it does not pick: it returns an owned experiment, so the
    /// copied records have to be copies. Source `pickExperiment` has the same
    /// two-map shape; a caller that does not need the input afterwards has no
    /// reason to pay for it.
    ///
    /// Unlike [`PeakPickerHiRes::pick_experiment`] this is **not atomic**: an
    /// error leaves the records before the failing one centroided and the rest
    /// as they were, and the returned report is lost. A caller that needs the
    /// input intact on failure uses [`PeakPickerHiRes::pick_experiment`].
    ///
    /// # Errors
    ///
    /// As [`PeakPickerHiRes::pick_experiment`].
    pub fn pick_experiment_in_place(
        &self,
        experiment: &mut MSExperiment,
    ) -> Result<PickedExperimentReport> {
        let limits = self.type_query_limits();
        let mut copies = self.start_experiment(experiment)?;
        let mut report = PickedExperimentReport {
            spectrum_boundaries: Vec::with_capacity(experiment.spectra.len()),
            chromatogram_boundaries: Vec::with_capacity(experiment.chromatograms.len()),
            omitted_spectrum_arrays: Vec::with_capacity(experiment.spectra.len()),
            omitted_chromatogram_arrays: Vec::with_capacity(experiment.chromatograms.len()),
        };
        for spectrum in &mut experiment.spectra {
            if !self.selects(spectrum, limits)? {
                report.spectrum_boundaries.push(None);
                report.omitted_spectrum_arrays.push(Vec::new());
                continue;
            }
            let picked = self.pick_spectrum_with_acquisition(spectrum, true, &mut copies)?;
            *spectrum = picked.spectrum;
            report.spectrum_boundaries.push(Some(picked.boundaries));
            report.omitted_spectrum_arrays.push(picked.omitted_arrays);
        }
        for chromatogram in &mut experiment.chromatograms {
            let picked =
                self.pick_chromatogram_with_acquisition(chromatogram, false, &mut copies)?;
            *chromatogram = picked.chromatogram;
            report.chromatogram_boundaries.push(picked.boundaries);
            report
                .omitted_chromatogram_arrays
                .push(picked.omitted_arrays);
        }
        Ok(report)
    }

    /// The spectrum-type query limits the experiment entry points use.
    fn type_query_limits(&self) -> SpectrumTypeQueryLimits {
        SpectrumTypeQueryLimits {
            max_points: self.max_points,
            ..SpectrumTypeQueryLimits::default()
        }
    }

    /// Validate the picker and the experiment and open its acquisition ledger.
    fn start_experiment(&self, input: &MSExperiment) -> Result<super::AcquisitionCopies> {
        self.validate()?;
        input.validate()?;
        let records = input
            .spectra
            .len()
            .checked_add(input.chromatograms.len())
            .ok_or_else(|| bad("the experiment record count overflows"))?;
        let mut copies = self.acquisition_ledger(records)?;
        copies.experiment(input)?;
        Ok(copies)
    }

    /// The acquisition-metadata copy ledger for an experiment of `records`
    /// records.
    ///
    /// The fixed part is the shared `AcquisitionCopies` default every other
    /// processing filter uses, read from that default rather than restated here
    /// so the two cannot drift;
    /// [`max_metadata_per_record`](Self::max_metadata_per_record) is added to
    /// both once per record, so the budget follows the input instead of putting
    /// a ceiling on how many records an experiment may have. The allowance is
    /// pooled rather than charged per record, as the shared ledger is: one
    /// record may spend another's share, and an input whose total acquisition
    /// metadata outweighs its own size is still refused.
    fn acquisition_ledger(&self, records: usize) -> Result<super::AcquisitionCopies> {
        let overflow = || bad("the acquisition metadata budget overflows for this record count");
        let allowance = records
            .checked_mul(self.max_metadata_per_record)
            .ok_or_else(overflow)?;
        let base = super::AcquisitionCopies::default();
        Ok(super::AcquisitionCopies {
            work: base.work.checked_add(allowance).ok_or_else(overflow)?,
            bytes: base.bytes.checked_add(allowance).ok_or_else(overflow)?,
        })
    }

    /// Whether the experiment entry points pick this spectrum, source
    /// `pickExperiment`'s automatic and manual mode selection.
    fn selects(&self, spectrum: &MSSpectrum, limits: SpectrumTypeQueryLimits) -> Result<bool> {
        if self.ms_levels.is_empty() {
            return Ok(spectrum.get_type_with_limits(true, limits)? != SpectrumType::Centroid);
        }
        if !self.ms_levels.contains(&spectrum.ms_level) {
            return Ok(false);
        }
        if spectrum.get_type_with_limits(true, limits)? == SpectrumType::Centroid
            && self.check_spectrum_type
        {
            return Err(bad(CENTROIDED_INPUT_MESSAGE));
        }
        Ok(true)
    }

    /// Replace a chromatogram with its picked form, only after picking succeeds.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerHiRes::pick_chromatogram`]; the input is unchanged then.
    pub fn filter_chromatogram(&self, input: &mut MSChromatogram) -> Result<()> {
        *input = self.pick_chromatogram(input)?.chromatogram;
        Ok(())
    }

    /// The source parameter defaults, as the `PeakPickerHiRes` constructor sets
    /// them: `signal_to_noise`, `spacing_difference_gap`, `spacing_difference`,
    /// `missing`, `ms_levels`, `report_FWHM`, `report_FWHM_unit`,
    /// `allow_missing_flank` and the `SignalToNoise:` subsection, with their
    /// value types, descriptions, restrictions and `advanced` tags, in
    /// declaration order.
    ///
    /// The boolean parameters are strings restricted to `true` and `false`, as
    /// in the source; a parameter file writer shows them as `bool`.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this picker in the layout of
    /// [`PeakPickerHiRes::defaults`], carrying the current values (source
    /// `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a value has no source parameter
    /// representation: `missing` or an MS level above `i32::MAX`, or a noise
    /// estimator value described at [`SignalToNoiseEstimatorMedian::to_param`].
    pub fn to_param(&self) -> Result<Param> {
        let missing = i32::try_from(self.missing)
            .map_err(|_| bad("'missing' does not fit the source 32-bit parameter"))?;
        let ms_levels = self
            .ms_levels
            .iter()
            .map(|&level| {
                i32::try_from(level)
                    .map_err(|_| bad("an MS level does not fit the source 32-bit parameter"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut levels = ParamEntry {
            name: "ms_levels".into(),
            description: "List of MS levels for which the peak picking is applied. If empty, auto mode is enabled, all peaks which aren't picked yet will get picked. Other scans are copied to the output without changes.".into(),
            value: ParamValue::IntegerList(ms_levels),
            ..ParamEntry::default()
        };
        levels.min_int = 1;
        let unit = self.report_fwhm.unwrap_or(self.inactive_fwhm_unit);
        let entries = vec![
            float_entry(
                "signal_to_noise",
                self.signal_to_noise,
                "Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)",
                false,
                Some(0.0),
                None,
            ),
            float_entry(
                "spacing_difference_gap",
                self.spacing_difference_gap,
                "The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.",
                true,
                Some(0.0),
                None,
            ),
            float_entry(
                "spacing_difference",
                self.spacing_difference,
                "Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.",
                true,
                Some(0.0),
                None,
            ),
            int_entry(
                "missing",
                i64::from(missing),
                "Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.",
                true,
                Some(0),
                None,
            ),
            levels,
            flag_entry(
                "report_FWHM",
                self.report_fwhm.is_some(),
                "Add metadata for FWHM (as floatDataArray named 'FWHM' or 'FWHM_ppm', depending on param 'report_FWHM_unit') for each picked peak.",
                false,
            ),
            string_entry(
                "report_FWHM_unit",
                match unit {
                    FwhmUnit::Absolute => "absolute",
                    FwhmUnit::Ppm => "relative",
                },
                "Unit of FWHM. Either absolute in the unit of input, e.g. 'm/z' for spectra, or relative as ppm (only sensible for spectra, not chromatograms).",
                false,
                &["relative", "absolute"],
            ),
            flag_entry(
                "allow_missing_flank",
                self.allow_missing_flank,
                "Allow peaks without flanking data points on both sides. This is useful for TimsTOF data where profile peaks may be missing the leading or trailing edge.",
                true,
            ),
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: vec![ParamNode {
                name: "SignalToNoise".into(),
                description: String::new(),
                entries: self.noise_estimator.param_entries()?,
                nodes: Vec::new(),
            }],
        })
    }

    /// Build a picker from parameters, as source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing parameters take their defaults, including the whole
    /// `SignalToNoise:` subsection. The native fields
    /// ([`check_spectrum_type`](Self::check_spectrum_type),
    /// [`ion_mobility_array`](Self::ion_mobility_array),
    /// [`compatibility`](Self::compatibility) and the limits) keep their
    /// defaults. Unknown parameters are ignored, as the source only warns; see
    /// [`PeakPickerHiRes::from_param_with_warnings`].
    ///
    /// `report_FWHM_unit` other than `absolute` means ppm, as the source's
    /// `!= "absolute"` test does. A `spacing_difference` or
    /// `spacing_difference_gap` of zero is kept as zero; the picker treats it
    /// as the source's infinity.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`PeakPickerHiRes::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerHiRes::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(PEAK_PICKER_HI_RES_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let unit = if param.value("report_FWHM_unit")?.as_str()? != "absolute" {
            FwhmUnit::Ppm
        } else {
            FwhmUnit::Absolute
        };
        let report = flag(param, "report_FWHM")?;
        let ms_levels = param
            .value("ms_levels")?
            .as_integer_list()?
            .iter()
            .map(|&level| u32::try_from(level).map_err(|_| bad("MS levels must not be negative")))
            .collect::<Result<Vec<_>>>()?;
        let defaults = Self::default();
        Ok(Self {
            signal_to_noise: float(param, "signal_to_noise")?,
            noise_estimator: SignalToNoiseEstimatorMedian::from_complete_param(
                &param.copy("SignalToNoise:", true)?,
            )?,
            spacing_difference: float(param, "spacing_difference")?,
            spacing_difference_gap: float(param, "spacing_difference_gap")?,
            missing: usize::try_from(integer(param, "missing")?)
                .map_err(|_| bad("'missing' must not be negative"))?,
            allow_missing_flank: flag(param, "allow_missing_flank")?,
            report_fwhm: report.then_some(unit),
            inactive_fwhm_unit: if report {
                defaults.inactive_fwhm_unit
            } else {
                unit
            },
            ms_levels,
            ..defaults
        })
    }
}

impl SpectrumFilter for PeakPickerHiRes {
    fn filter_spectrum(&self, input: &mut MSSpectrum) -> Result<()> {
        *input = self.pick_spectrum(input)?.spectrum;
        Ok(())
    }
    fn filter_experiment(&self, input: &mut MSExperiment) -> Result<()> {
        *input = self.pick_experiment(input)?.experiment;
        Ok(())
    }
}

fn fwhm_name(unit: FwhmUnit) -> &'static str {
    match unit {
        FwhmUnit::Absolute => "FWHM",
        FwhmUnit::Ppm => "FWHM_ppm",
    }
}

/// Infer whether a spectrum holds profile or centroid data from up to five
/// high-intensity peak shoulders.
///
/// The source `PeakTypeEstimator` heuristic: fewer than five samples give
/// [`SpectrumType::Unknown`], and no positive evidence gives
/// [`SpectrumType::Centroid`]. Unlike [`MSSpectrum::get_type`], this ignores the
/// stored type and processing history and always inspects the data.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] for an invalid spectrum, more than one
/// million samples, non-finite values, negative intensities or duplicate
/// positions, and [`Error::UnsortedData`] for decreasing positions. The source
/// heuristic classifies such data without checking; this entry point keeps the
/// picker's stricter input contract.
pub fn estimate_spectrum_type(input: &MSSpectrum) -> Result<SpectrumType> {
    estimate_spectrum_type_with_limit(input, 1_000_000)
}
fn estimate_spectrum_type_with_limit(
    input: &MSSpectrum,
    max_points: usize,
) -> Result<SpectrumType> {
    if input.len() > max_points {
        return Err(bad("spectrum type estimation exceeds point limit"));
    }
    input.validate()?;
    let x: Vec<_> = input.peaks.iter().map(|p| p.mz).collect();
    let mut y: Vec<_> = input.peaks.iter().map(|p| f64::from(p.intensity)).collect();
    validate_points(
        x.len(),
        &|i| x[i],
        &|i| y[i],
        max_points,
        &PickingCompatibility::default(),
    )?;
    Ok(crate::kernel::spectrum_type::estimate(&x, &mut y))
}

#[cfg(test)]
mod acquisition_ledger_tests {
    use super::PeakPickerHiRes;
    use crate::Error;
    use crate::processing::AcquisitionCopies;

    /// The ledger the picker opens for an experiment of `records` records, as
    /// its remaining work units and bytes.
    ///
    /// `AcquisitionCopies` is a shared type in `src/processing.rs` and does not
    /// implement `Debug`, so the ledger is unwrapped by hand rather than with
    /// `Result::unwrap`.
    fn ledger(picker: &PeakPickerHiRes, records: usize) -> (usize, usize) {
        match picker.acquisition_ledger(records) {
            Ok(opened) => (opened.work, opened.bytes),
            Err(error) => panic!("{records} records were refused a ledger: {error:?}"),
        }
    }

    #[test]
    fn the_acquisition_ledger_follows_the_record_count() {
        // The shared ledger's fixed part is a ceiling on how many records an
        // experiment may have, which source `pickExperiment` has no counterpart
        // for: an ordinary vendor-converted run spends about 8 KiB of it per
        // spectrum, so it stops at roughly 34 000 spectra and refused the
        // 40 856 of the 2.3 GB Q Exactive benchmark run outright.
        // `max_metadata_per_record` is added to both dimensions once per input
        // record, so the budget follows the input instead of capping it.
        //
        // This reads the ledger the picker opens rather than building an
        // experiment whose metadata exhausts 256 MiB: the meter charges a
        // record within about a factor of two of what that record actually
        // costs in memory, so crossing the fixed part end to end costs hundreds
        // of mebibytes of test process. The end-to-end behaviour of the
        // allowance is pinned by `tests/peak_picking_experiment.rs`, and the
        // real 40 856-spectrum run is evidence in
        // `docs/PEAK_PICKING_SUPPORT.md`.
        let base = AcquisitionCopies::default();
        let picker = PeakPickerHiRes::default();
        assert_eq!(picker.max_metadata_per_record, 64 * 1024);
        assert_eq!(ledger(&picker, 0), (base.work, base.bytes));
        for records in [1_usize, 2, 34_257, 40_856, 1 << 20] {
            let allowance = records * picker.max_metadata_per_record;
            assert_eq!(
                ledger(&picker, records),
                (base.work + allowance, base.bytes + allowance),
                "{records} records"
            );
        }
        // Strictly increasing in the record count, in both dimensions, so no
        // record count is a ceiling.
        assert!(ledger(&picker, 40_857) > ledger(&picker, 40_856));
        // Zero pins the fixed part exactly, whatever the record count: the
        // behaviour before the field existed.
        let pinned = PeakPickerHiRes {
            max_metadata_per_record: 0,
            ..Default::default()
        };
        for records in [0_usize, 1, 40_856, usize::MAX] {
            assert_eq!(ledger(&pinned, records), (base.work, base.bytes));
        }
    }

    #[test]
    fn an_overflowing_metadata_allowance_is_refused() {
        // Both the multiplication by the record count and the addition to the
        // fixed part are checked, so an absurd allowance is refused instead of
        // wrapping into a budget below the fixed one.
        for allowance in [usize::MAX, usize::MAX / 2, usize::MAX / 4] {
            let picker = PeakPickerHiRes {
                max_metadata_per_record: allowance,
                ..Default::default()
            };
            match picker.acquisition_ledger(8) {
                Ok(_) => panic!("an allowance of {allowance} per record was admitted"),
                Err(Error::InvalidValue(text)) => assert!(text.contains("overflows"), "{text}"),
                Err(other) => panic!("{other:?}"),
            }
        }
        // One record of the largest possible allowance overflows only in the
        // addition to the fixed part, which is checked as well.
        let picker = PeakPickerHiRes {
            max_metadata_per_record: usize::MAX,
            ..Default::default()
        };
        assert!(picker.acquisition_ledger(1).is_err());
        // A zero allowance cannot overflow, however many records.
        let pinned = PeakPickerHiRes {
            max_metadata_per_record: 0,
            ..Default::default()
        };
        assert!(pinned.acquisition_ledger(usize::MAX).is_ok());
    }
}
