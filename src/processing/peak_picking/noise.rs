// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Histogram-median signal-to-noise estimation.
//!
//! Ports the header-only template
//! `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h`.
//! Its base `SignalToNoiseEstimator.h` is
//! [`crate::processing::noise_estimation`]. `docs/SIGNAL_TO_NOISE_SUPPORT.md`
//! holds the API mapping of every public source member, the preserved
//! conventions, the native differences and the evidence.
//!
//! The C++ class is a `DefaultParamHandler` and `ProgressLogger` that computes
//! its estimates in `init` and serves them through `getSignalToNoise(index)`.
//! This port is stateless: [`SignalToNoiseEstimatorMedian::estimate`] and its
//! variants return all estimates, both percentages and the warnings at once,
//! and the parameter contract maps through `defaults`, `from_param`,
//! `set_parameters` and `to_param` on the same type.
//!
//! # Native safety profile and source behaviour
//!
//! [`PickingCompatibility::default`] is the native safety profile: it refuses
//! inputs and parameter values on which the source computes non-finite or
//! degenerate results, and it bins out-of-range intensities by the documented
//! intent. [`PickingCompatibility::source`], which every TOPP tool path uses,
//! computes what the Linux x86-64 Release build computes wherever the source is
//! defined, and refuses exactly where the source is undefined.

use super::{PickingCompatibility, SignalPoint, bad};
use crate::concept::progress_logger::ProgressLogger;
use crate::kernel::{ChromatogramPeak, MSChromatogram, MSSpectrum, Peak1D};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::processing::noise_estimation::{GaussianEstimate, SignalToNoiseEstimator, x86};
use crate::{Error, Result};

/// How the upper end of the intensity histogram is chosen.
///
/// Source `SignalToNoiseEstimatorMedian::IntensityThresholdCalculation` ("method
/// to use for estimating the maximal intensity that is used for histogram
/// calculation") with the parameter `auto_mode`: `MANUAL = -1`,
/// `AUTOMAXBYSTDEV = 0` (the default) and `AUTOMAXBYPERCENT = 1`. Intensities
/// at or above the upper end go into the last histogram bin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NoiseHistogramRange {
    /// `auto_mode = 0`: mean plus `factor` times the population standard
    /// deviation of all intensities (`auto_max_stdev_factor`, source range
    /// `0..=999`; a NaN passes the source restriction and is accepted under
    /// [`NoiseCompatibility::source_value_domain`]).
    StandardDeviation {
        /// Multiplier of the standard deviation.
        factor: f64,
    },
    /// `auto_mode = -1`: a fixed upper end (`max_intensity`). The source throws
    /// `Exception::InvalidValue` when estimation runs with a value `<= 0`; this
    /// returns [`Error::InvalidValue`] with the source text at the same point.
    Manual {
        /// Upper end of the histogram. The source parameter is an integer; any
        /// positive finite value is accepted here.
        max_intensity: f64,
    },
    /// `auto_mode = 1`: the `auto_max_percentile`-th percentile of a 100-bin
    /// pre-histogram, computed exactly as the source computes it on the inputs
    /// where the source is defined.
    ///
    /// The source (`SignalToNoiseEstimatorMedian.h:191-233`) finds the
    /// **minimum** intensity `m` (the comparator `a > b` at `:208` reverses
    /// `std::max_element`), sets `bin_size = m / 100` in `f32` (`:211`),
    /// increments the pre-histogram at `(int)((I - 1) / bin_size)`, with
    /// `I - 1` in `f32`, without a bounds check (`:216`), and walks at most one
    /// bin per container element (`:225-230`), so the walk ends after `n`
    /// bins even when the cumulative count has not reached the target. The
    /// upper end is `(i + 0.5) * bin_size`, where `i` is the last bin visited
    /// (`-1` when none is), which a target of zero makes negative.
    ///
    /// The source is defined exactly on non-empty containers in which every
    /// intensity's quotient `q = (I - 1) / bin_size` satisfies `-1 < q < 100`
    /// and whose `int` counters do not overflow: then every index lies in
    /// `[0, 99]`, the conversion is in range, and the walk stops by bin 99
    /// because the 100 bins hold all `n` elements. Roughly, the minimum must
    /// exceed `100 / 101` and every intensity must lie below `m + 1`. The
    /// counters only matter beyond `i32::MAX` points: every pre-histogram bin
    /// must hold at most `i32::MAX` points (`:216`), and the walk's running
    /// count must stay within `int` (`:228`). The walk's target
    /// `p * n / 100` (`:220`) does not fit `int` from `2^31` on; the Release
    /// build's 32-bit `cvttsd2si` then returns `INT_MIN`, which skips the walk,
    /// and this port reproduces that measured outcome. A target of zero or
    /// `INT_MIN` makes the upper end negative, so the source warns and returns
    /// before its main loop, whatever the number of points; any other upper end
    /// reaches the main loop, whose counters take at most `i32::MAX` points.
    /// Anywhere else estimation returns [`Error::Unsupported`] naming the source
    /// line that becomes undefined: `:209` dereferences `end()` of an empty
    /// container, `:216` writes outside the pre-histogram or overflows a bin
    /// count, `:228` overflows the running count, and the main loop overflows
    /// its window count (`:365`). For the out-of-range writes there is no
    /// answer to reproduce: of nine such inputs run three times
    /// each on the Release build (`../oracle/sne-completion/probes/`), seven
    /// end with `SIGSEGV` or `SIGABRT`, and two, a write one bin past the end
    /// and a negative minimum, return whatever the neighbouring memory then
    /// yields; the product SDK ended with `SIGBUS` or `SIGSEGV` (CPP-256).
    /// `docs/SIGNAL_TO_NOISE_SUPPORT.md` writes the domain derivation out.
    Percentile {
        /// The percentile, `0..=100` in the source.
        percentile: f64,
    },
}
impl Default for NoiseHistogramRange {
    fn default() -> Self {
        Self::StandardDeviation { factor: 3.0 }
    }
}

/// The three histogram-range parameter values, whether or not the selected
/// [`NoiseHistogramRange`] uses them.
///
/// The source object keeps every parameter it was given, including those the
/// selected `auto_mode` ignores, and returns them from `getParameters`. This
/// record keeps the ignored ones so that
/// [`SignalToNoiseEstimatorMedian::to_param`] reproduces the parameters
/// [`SignalToNoiseEstimatorMedian::from_param`] received. The value the active
/// mode uses is taken from [`SignalToNoiseEstimatorMedian::histogram_range`]
/// instead of from here.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NoiseRangeParameters {
    /// `max_intensity`, source default `-1`.
    pub max_intensity: i32,
    /// `auto_max_stdev_factor`, source default `3.0`.
    pub auto_max_stdev_factor: f64,
    /// `auto_max_percentile`, source default `95`.
    pub auto_max_percentile: i32,
}
impl Default for NoiseRangeParameters {
    fn default() -> Self {
        Self {
            max_intensity: -1,
            auto_max_stdev_factor: 3.0,
            auto_max_percentile: 95,
        }
    }
}

/// How an intensity's histogram quotient becomes a bin index.
///
/// The source converts before it clamps
/// (`std::max(std::min<int>((int)(I / bin_size), bin_count - 1), 0)` at
/// `SignalToNoiseEstimatorMedian.h:297` and `:308`), so a quotient outside the
/// `int` range, or a NaN, is undefined behaviour (CPP-257). Such quotients arise
/// with legal parameters: a manual `max_intensity` of `1` and an intensity above
/// `2^31`, a standard-deviation `factor` of `0` with a large `bin_count`, or
/// negative intensities.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BinIndexConversion {
    /// Clamp the quotient to `[0, bin_count - 1]` before truncating: an
    /// intensity at or above the histogram's upper end lands in the last bin,
    /// as the `max_intensity` parameter description promises, and a negative
    /// quotient in the first. A NaN quotient, which only the source value
    /// domain admits, lands in the last bin. For every non-NaN quotient this is
    /// also what the arm64 build computes, whose `fcvtzs` saturates (the P1
    /// verifier's `huge_int` case); for NaN, `fcvtzs` gives `0`. This is the
    /// default.
    #[default]
    ClampBeforeTruncation,
    /// The Linux x86-64 Release build: the 32-bit `cvttsd2si` both
    /// `libOpenMS.so` instantiations emit returns `INT_MIN` for every quotient
    /// whose truncation does not fit and for NaN, which the clamp then sends to
    /// bin **0**. Measured bit for bit by the `cpp257_*` oracle cases, including
    /// `PeakPickerHiRes::pick`.
    X86_64Release,
}

/// The source behaviours of the median estimator that the native safety profile
/// refuses or computes differently.
///
/// Part of [`PickingCompatibility`]; [`NoiseCompatibility::source`] selects all
/// of them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NoiseCompatibility {
    /// Compute with the values the source computes with instead of refusing
    /// them: the parameter values its restrictions accept (a NaN or infinite
    /// `win_len`, a NaN `auto_max_stdev_factor`, a zero, negative, infinite or
    /// NaN `noise_for_empty_window`), non-finite positions and intensities, and
    /// the infinite or NaN window bounds, histogram ranges and ratios IEEE
    /// arithmetic then produces (a tiny `noise_for_empty_window` stores an
    /// infinite ratio). NaN results carry the bits the Linux x86-64 Release
    /// build produces.
    pub source_value_domain: bool,
    /// For an empty input return NaN percentages, from the source's
    /// `0 * 100 / 0` (`SignalToNoiseEstimatorMedian.h:373-374`), and in the
    /// standard-deviation range a NaN
    /// `max_intensity`, from `0 / 0` in `SignalToNoiseEstimator::estimate_`
    /// (`SignalToNoiseEstimator.h:127`); both are the x86-64 default NaN
    /// (`0xfff8000000000000`). The native profile returns zero for all three.
    pub nan_for_empty_input: bool,
    /// How an out-of-range histogram quotient is binned.
    pub bin_index: BinIndexConversion,
}
impl NoiseCompatibility {
    /// Every source behaviour, as the Linux x86-64 Release build computes it.
    pub const fn source() -> Self {
        Self {
            source_value_domain: true,
            nan_for_empty_input: true,
            bin_index: BinIndexConversion::X86_64Release,
        }
    }
}

/// A data point the estimator reads: the source's `PeakType`
/// (`Container::const_iterator::value_type`) reduced to what `computeSTN_`
/// uses.
pub trait NoisePoint {
    /// The position (`getPos()`): m/z for a spectrum, retention time in seconds
    /// for a chromatogram.
    fn position(&self) -> f64;
    /// The stored intensity (`getIntensity()`, a `float`).
    fn intensity(&self) -> f32;
}
impl NoisePoint for Peak1D {
    fn position(&self) -> f64 {
        self.mz
    }
    fn intensity(&self) -> f32 {
        self.intensity
    }
}
impl NoisePoint for ChromatogramPeak {
    fn position(&self) -> f64 {
        self.rt
    }
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

/// Estimates the signal-to-noise ratio of every data point from the median of
/// a sliding-window intensity histogram.
///
/// For each data point, the points within half of
/// [`window_length`](Self::window_length) on either side (both ends inclusive)
/// form its window. The noise is the median intensity of that window, read
/// from a histogram of [`bin_count`](Self::bin_count) bins over
/// `[0, upper end]`, which determines the level of error and the runtime: the
/// median bin is the first whose cumulative count reaches `(count + 1) / 2`, and
/// the noise is linearly interpolated inside that bin, assuming its points are
/// uniformly distributed, then floored at one. Interpolating keeps the estimate
/// close to continuous as the data changes, where the bin centre alone could
/// take only `bin_count` values. Bins are at least one intensity unit wide. A
/// window with fewer than
/// [`min_required_elements`](Self::min_required_elements) points is *sparse*
/// and uses [`noise_for_empty_window`](Self::noise_for_empty_window). The
/// signal-to-noise ratio is the point's intensity divided by its noise. The
/// upper end is determined automatically by one of two methods or set by the
/// caller ([`NoiseHistogramRange`]).
///
/// # Notes
///
/// If the upper end is too low and the median falls into the last bin, the
/// estimate is unreliable: increase `max_intensity`, and possibly `bin_count`.
/// The source warns through `OPENMS_LOG_WARN` when more than 20 percent of
/// windows are sparse or more than 1 percent of medians fall into the last bin,
/// unless [`write_log_messages`](Self::write_log_messages) is off; it always
/// warns when an automatic upper end comes out negative. These lines are
/// returned in [`NoiseEstimates::log`], and both percentages in
/// [`NoiseEstimates`].
///
/// Changing a parameter invalidates the source's stored estimates, which it
/// then recomputes on the next request; this port recomputes on every call and
/// has no cache to invalidate. The source's copy constructor and assignment
/// clear the copy's estimates for the same reason; a [`Clone`] here has none.
#[derive(Clone, Debug, PartialEq)]
pub struct SignalToNoiseEstimatorMedian {
    /// Upper end of the intensity histogram (`auto_mode` and the parameter it
    /// selects).
    pub histogram_range: NoiseHistogramRange,
    /// Values of the range parameters the selected mode does not use; see
    /// [`NoiseRangeParameters`].
    pub range_parameters: NoiseRangeParameters,
    /// Full window width `win_len` ("window length in Thomson") in coordinate
    /// units (Th for spectra, seconds for chromatograms). At least one, or NaN
    /// or infinite under [`NoiseCompatibility::source_value_domain`]; source
    /// default `200`. An infinite length makes every window the whole input.
    pub window_length: f64,
    /// Number of histogram bins `bin_count`, at least three; source default
    /// `30`.
    pub bin_count: usize,
    /// Minimum points in a window, `min_required_elements`, at least one;
    /// source default `10`.
    pub min_required_elements: usize,
    /// Noise used for sparse windows, `noise_for_empty_window`; source default
    /// `1e20` ("use a very high value if you want to get a low S/N result").
    /// The native profile requires a finite positive value; the source accepts
    /// any.
    pub noise_for_empty_window: f64,
    /// Source `write_log_messages`: whether the sparse-window and rightmost-bin
    /// warnings are returned in [`NoiseEstimates::log`]. The negative-range
    /// warning is returned regardless, as in the source.
    pub write_log_messages: bool,
    /// Native ceiling on the number of input points.
    pub max_points: usize,
    /// Native ceiling on [`bin_count`](Self::bin_count).
    pub max_bins: usize,
    /// Native ceiling on histogram insertions, removals and median-bin visits.
    pub max_work: usize,
}
impl Default for SignalToNoiseEstimatorMedian {
    fn default() -> Self {
        Self {
            histogram_range: Default::default(),
            range_parameters: Default::default(),
            window_length: 200.0,
            bin_count: 30,
            min_required_elements: 10,
            noise_for_empty_window: 1e20,
            write_log_messages: true,
            max_points: 1_000_000,
            max_bins: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

/// Result of [`SignalToNoiseEstimatorMedian::estimate`].
#[derive(Clone, Debug, PartialEq)]
pub struct NoiseEstimates {
    /// Signal-to-noise ratio per input point, source `getSignalToNoise(i)`
    /// (the member `stn_estimates_`).
    pub signal_to_noise: Vec<f64>,
    /// Noise per input point (native). `f64::INFINITY` for every point when the
    /// source returns early; see [`SignalToNoiseEstimatorMedian::estimate`].
    pub noise: Vec<f64>,
    /// The histogram upper end that was used (the source member
    /// `max_intensity_` after estimation).
    pub max_intensity: f64,
    /// Percentage of sparse windows, source `getSparseWindowPercent` ("how many
    /// percent of the windows were sparse").
    pub sparse_window_percent: f64,
    /// Percentage of medians found in the last bin, source
    /// `getHistogramRightmostPercent`.
    pub histogram_rightmost_percent: f64,
    /// The lines the source writes to `OPENMS_LOG_WARN` during this
    /// estimation, without their line ends, in order. The source's log stream
    /// then suppresses repeated lines and reports them as
    /// `<line> occurred N times`, which `crate::concept::log_stream` ports.
    pub log: Vec<String>,
}

/// Parameter handler name of the source class, used in parameter diagnostics.
pub const SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME: &str = "SignalToNoiseEstimatorMedian";

/// The label of the estimator's progress report, source
/// `startProgress(0, c.size(), "noise estimation of data")`.
pub const NOISE_PROGRESS_LABEL: &str = "noise estimation of data";

const MAX_INTENSITY_DESCRIPTION: &str = "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin. If you choose 'max_intensity' too small, the noise estimate might be too small as well.  If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime). In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.";

/// The largest count the source's `int` counters hold. Beyond this many
/// points, `estimate_` (`size`, `SignalToNoiseEstimator.h:123`) and the main
/// loop of `computeSTN_` (`window_count` at `SignalToNoiseEstimatorMedian.h:365`)
/// overflow; AUTOMAXBYPERCENT's
/// counters (`:216`, `:228`) overflow only when a pre-histogram bin or the
/// walk's running count exceeds it, and a negative range returns before the
/// main loop. Each refusal is placed where its counter overflows.
const SOURCE_MAX_POINTS: usize = i32::MAX as usize;

fn entry(name: &str, value: ParamValue, description: &str, advanced: bool) -> ParamEntry {
    let mut entry = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        entry.tags.insert("advanced".into());
    }
    entry
}

pub(super) fn int_entry(
    name: &str,
    value: i64,
    description: &str,
    advanced: bool,
    min: Option<i32>,
    max: Option<i32>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Integer(value), description, advanced);
    if let Some(min) = min {
        e.min_int = min;
    }
    if let Some(max) = max {
        e.max_int = max;
    }
    e
}

pub(super) fn float_entry(
    name: &str,
    value: f64,
    description: &str,
    advanced: bool,
    min: Option<f64>,
    max: Option<f64>,
) -> ParamEntry {
    let mut e = entry(name, ParamValue::Float(value), description, advanced);
    if let Some(min) = min {
        e.min_float = min;
    }
    if let Some(max) = max {
        e.max_float = max;
    }
    e
}

pub(super) fn flag_entry(name: &str, value: bool, description: &str, advanced: bool) -> ParamEntry {
    string_entry(
        name,
        if value { "true" } else { "false" },
        description,
        advanced,
        &["true", "false"],
    )
}

pub(super) fn string_entry(
    name: &str,
    value: &str,
    description: &str,
    advanced: bool,
    valid: &[&str],
) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(value.into()),
        description,
        advanced,
    );
    e.valid_strings = valid.iter().map(|v| (*v).to_owned()).collect();
    e
}

pub(super) fn integer(param: &Param, key: &str) -> Result<i64> {
    param.value(key)?.to_i64()
}

pub(super) fn float(param: &Param, key: &str) -> Result<f64> {
    match param.value(key)? {
        ParamValue::Float(value) => Ok(*value),
        _ => Err(bad(&format!(
            "parameter '{key}' must be a floating-point value"
        ))),
    }
}

pub(super) fn flag(param: &Param, key: &str) -> Result<bool> {
    param.value(key)?.to_bool()
}

fn to_i32(value: usize, key: &str) -> Result<i64> {
    i32::try_from(value)
        .map(i64::from)
        .map_err(|_| bad(&format!("'{key}' does not fit the source 32-bit parameter")))
}

fn integral_i32(value: f64, key: &str) -> Result<i32> {
    if value.fract() == 0.0 && value >= f64::from(i32::MIN) && value <= f64::from(i32::MAX) {
        Ok(value as i32)
    } else {
        Err(bad(&format!(
            "'{key}' = {value} is not representable as the source integer parameter"
        )))
    }
}

/// `Exception::InvalidValue(..., message, StringUtils::toStr(value))::what()`.
fn source_invalid_value(value: f64, message: &str) -> Error {
    Error::InvalidValue(format!(
        "the value '{}' was used but is not valid; {message}",
        crate::param::value::format_float(value, true)
    ))
}

/// `std::ostream << double` at the stream's default precision 6, which is C
/// `printf("%g")`: six significant digits, fixed notation when the decimal
/// exponent `X` after rounding satisfies `-4 <= X < 6` and scientific notation
/// with a signed two-digit exponent otherwise, trailing fraction zeros and a
/// bare point removed. Exact decimal ties round half to even, as glibc does.
///
/// This is the rule `format::file_info::text_format::ostream_g` implements at
/// precision 6. The processing module may not depend on the format module
/// (`docs/module-cycles.json`: that edge would close a cycle), so the warning
/// lines restate it here. Only finite values reach it: the negative-range
/// warning prints a range below zero, which is never NaN or `-inf` (a sum that
/// reaches `-inf` makes the variance, and so the range, NaN), and the other two
/// print a percentage above 1 or 20. Non-finite values still get glibc's
/// spellings.
fn stream_double(value: f64) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .into();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.into();
    }
    fn strip(text: &str) -> &str {
        if text.contains('.') {
            text.trim_end_matches('0').trim_end_matches('.')
        } else {
            text
        }
    }
    // Rust's exact formatting rounds half to even, so this exponent is the one
    // after rounding to six significant digits.
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i32 = exponent.parse().unwrap_or(0);
    if (-4..6).contains(&exponent) {
        // 0 <= 5 - exponent <= 9.
        let fraction = usize::try_from(5 - exponent).unwrap_or(0);
        strip(&format!("{value:.fraction$}")).to_owned()
    } else {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", strip(mantissa), exponent.unsigned_abs())
    }
}

const MANUAL_MESSAGE: &str = "auto_mode is on MANUAL! max_intensity is <=0. Needs to be positive! Use setMaxIntensity(<value>) or enable auto_mode!";
const PERCENTILE_MESSAGE: &str = "auto_mode is on AUTOMAXBYPERCENT! auto_max_percentile is not in [0,100]. Use setAutoMaxPercentile(<value>) to change it!";

impl SignalToNoiseEstimatorMedian {
    /// The source parameter defaults, as `SignalToNoiseEstimatorMedian()` sets
    /// them: names, values, value types, descriptions, restrictions and the
    /// `advanced` tag, in declaration order.
    ///
    /// Source `getDefaults()`, the public member `defaults_` (re-exported by
    /// `using` at `SignalToNoiseEstimatorMedian.h:71`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    pub(super) fn param_entries(&self) -> Result<Vec<ParamEntry>> {
        let r = &self.range_parameters;
        let (auto_mode, max_intensity, factor, percentile) = match self.histogram_range {
            NoiseHistogramRange::Manual { max_intensity } => (
                -1,
                integral_i32(max_intensity, "max_intensity")?,
                r.auto_max_stdev_factor,
                r.auto_max_percentile,
            ),
            NoiseHistogramRange::StandardDeviation { factor } => {
                (0, r.max_intensity, factor, r.auto_max_percentile)
            }
            NoiseHistogramRange::Percentile { percentile } => (
                1,
                r.max_intensity,
                r.auto_max_stdev_factor,
                integral_i32(percentile, "auto_max_percentile")?,
            ),
        };
        Ok(vec![
            int_entry(
                "max_intensity",
                i64::from(max_intensity),
                MAX_INTENSITY_DESCRIPTION,
                true,
                Some(-1),
                None,
            ),
            float_entry(
                "auto_max_stdev_factor",
                factor,
                "parameter for 'max_intensity' estimation (if 'auto_mode' == 0): mean + 'auto_max_stdev_factor' * stdev",
                true,
                Some(0.0),
                Some(999.0),
            ),
            int_entry(
                "auto_max_percentile",
                i64::from(percentile),
                "parameter for 'max_intensity' estimation (if 'auto_mode' == 1): auto_max_percentile th percentile",
                true,
                Some(0),
                Some(100),
            ),
            int_entry(
                "auto_mode",
                auto_mode,
                "method to use to determine maximal intensity: -1 --> use 'max_intensity'; 0 --> 'auto_max_stdev_factor' method (default); 1 --> 'auto_max_percentile' method",
                true,
                Some(-1),
                Some(1),
            ),
            float_entry(
                "win_len",
                self.window_length,
                "window length in Thomson",
                false,
                Some(1.0),
                None,
            ),
            int_entry(
                "bin_count",
                to_i32(self.bin_count, "bin_count")?,
                "number of bins for intensity values",
                false,
                Some(3),
                None,
            ),
            int_entry(
                "min_required_elements",
                to_i32(self.min_required_elements, "min_required_elements")?,
                "minimum number of elements required in a window (otherwise it is considered sparse)",
                false,
                Some(1),
                None,
            ),
            float_entry(
                "noise_for_empty_window",
                self.noise_for_empty_window,
                "noise value used for sparse windows",
                true,
                None,
                None,
            ),
            string_entry(
                "write_log_messages",
                if self.write_log_messages {
                    "true"
                } else {
                    "false"
                },
                "Write out log messages in case of sparse windows or median in rightmost histogram bin",
                false,
                &["true", "false"],
            ),
        ])
    }

    /// The parameters of this estimator in the source layout of
    /// [`SignalToNoiseEstimatorMedian::defaults`], carrying the current values.
    ///
    /// The value side of source `getParameters()`: here the typed fields are
    /// the state and the tree is derived from them, so its names, types,
    /// descriptions, tags and restrictions are always the defaults'. The
    /// source member `param_` (public through the `using` at
    /// `SignalToNoiseEstimatorMedian.h:72`) is instead the caller's tree merged
    /// with the defaults, so it keeps what the caller's entries carried: a
    /// `Param::setValue(key, value)` without tags, for example, leaves
    /// `param_` without the `advanced` tag. That tree, exactly, is
    /// [`DefaultParamHandler::parameters`] of the handler
    /// [`SignalToNoiseEstimatorMedian::from_param_with_handler`] returns. The
    /// values of both trees always agree.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a value has no source parameter
    /// representation: a manual `max_intensity` or a percentile that is not an
    /// integer in the 32-bit range, or a `bin_count` or `min_required_elements`
    /// above `i32::MAX`.
    pub fn to_param(&self) -> Result<Param> {
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: self.param_entries()?,
            nodes: Vec::new(),
        })
    }

    /// Build an estimator from parameters, as source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing parameters take their defaults; the native resource ceilings keep
    /// their defaults. Unknown parameters are ignored, as the source only warns
    /// about them; use [`SignalToNoiseEstimatorMedian::from_param_with_warnings`]
    /// to see those warnings.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source `Exception::InvalidParameter`),
    /// or when `bin_count` or `min_required_elements` does not fit `usize`.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`SignalToNoiseEstimatorMedian::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let (estimator, _, warnings) = Self::from_param_with_handler(parameters)?;
        Ok((estimator, warnings))
    }

    /// As [`SignalToNoiseEstimatorMedian::from_param_with_warnings`], also
    /// returning the parameter handler a source object holds after
    /// `setParameters(parameters)`.
    ///
    /// The handler carries the source's public members `defaults_`
    /// ([`DefaultParamHandler::defaults`], equal to
    /// [`SignalToNoiseEstimatorMedian::defaults`]) and `param_`
    /// ([`DefaultParamHandler::parameters`], the given tree merged with the
    /// defaults, as `getParameters()` returns it). Both are public in the
    /// source through the `using` declarations at
    /// `SignalToNoiseEstimatorMedian.h:71-72`; as there, changing the tree
    /// afterwards does not change the estimator, whose members only
    /// `updateMembers_`, i.e. this constructor, derives.
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::from_param`].
    pub fn from_param_with_handler(
        parameters: &Param,
    ) -> Result<(Self, DefaultParamHandler, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        let (estimator, warnings) =
            handler.set_parameters_with(parameters, Self::from_complete_param)?;
        Ok((estimator, handler, warnings))
    }

    /// Replace the parameter values of this estimator, as source
    /// `setParameters` on an existing object, and return the unknown-parameter
    /// warnings.
    ///
    /// The native resource ceilings are kept. On error the estimator is
    /// unchanged.
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::from_param`].
    pub fn set_parameters(&mut self, parameters: &Param) -> Result<Vec<String>> {
        let (updated, warnings) = Self::from_param_with_warnings(parameters)?;
        *self = Self {
            max_points: self.max_points,
            max_bins: self.max_bins,
            max_work: self.max_work,
            ..updated
        };
        Ok(warnings)
    }

    /// Typed members from a parameter tree that already holds every default.
    pub(super) fn from_complete_param(param: &Param) -> Result<Self> {
        let range_parameters = NoiseRangeParameters {
            max_intensity: param.value("max_intensity")?.to_i32()?,
            auto_max_stdev_factor: float(param, "auto_max_stdev_factor")?,
            auto_max_percentile: param.value("auto_max_percentile")?.to_i32()?,
        };
        // Source computeSTN_: 0 is AUTOMAXBYSTDEV, 1 AUTOMAXBYPERCENT and every
        // other value MANUAL; the restriction `-1:1` leaves only -1 for the last.
        let histogram_range = match integer(param, "auto_mode")? {
            0 => NoiseHistogramRange::StandardDeviation {
                factor: range_parameters.auto_max_stdev_factor,
            },
            1 => NoiseHistogramRange::Percentile {
                percentile: f64::from(range_parameters.auto_max_percentile),
            },
            _ => NoiseHistogramRange::Manual {
                max_intensity: f64::from(range_parameters.max_intensity),
            },
        };
        let count = |key: &str| -> Result<usize> {
            usize::try_from(integer(param, key)?)
                .map_err(|_| bad(&format!("parameter '{key}' must not be negative")))
        };
        Ok(Self {
            histogram_range,
            range_parameters,
            window_length: float(param, "win_len")?,
            bin_count: count("bin_count")?,
            min_required_elements: count("min_required_elements")?,
            noise_for_empty_window: float(param, "noise_for_empty_window")?,
            write_log_messages: flag(param, "write_log_messages")?,
            ..Self::default()
        })
    }

    /// The option checks estimation runs first, in two groups: values the
    /// source's parameter restrictions refuse (or that its members cannot hold)
    /// and native refusals, then the two source exceptions of `computeSTN_`.
    fn validate(&self, compatibility: &NoiseCompatibility) -> Result<()> {
        let source = compatibility.source_value_domain;
        // `win_len >= 1` (Param.cpp skips the unset upper bound, so +inf and,
        // since `NaN < 1` is false, NaN pass the restriction).
        let window_ok = if source {
            self.window_length >= 1.0 || self.window_length.is_nan()
        } else {
            self.window_length.is_finite() && self.window_length >= 1.0
        };
        let empty_noise_ok = source
            || (self.noise_for_empty_window.is_finite() && self.noise_for_empty_window > 0.0);
        if !window_ok
            || self.bin_count < 3
            || self.bin_count > self.max_bins
            || self.bin_count > SOURCE_MAX_POINTS
            || self.min_required_elements == 0
            || !empty_noise_ok
            || self.max_points == 0
            || self.max_work == 0
        {
            return Err(bad("invalid median-noise options or resource limits"));
        }
        match self.histogram_range {
            NoiseHistogramRange::StandardDeviation { factor } => {
                let in_range = (0.0..=999.0).contains(&factor);
                if in_range || (source && factor.is_nan()) {
                    Ok(())
                } else {
                    Err(bad("noise standard-deviation factor must be in 0..=999"))
                }
            }
            // SignalToNoiseEstimatorMedian.h:236-244. The source member is an
            // integer parameter, so it is never NaN or infinite.
            NoiseHistogramRange::Manual { max_intensity }
                if !max_intensity.is_finite() || max_intensity <= 0.0 =>
            {
                Err(source_invalid_value(max_intensity, MANUAL_MESSAGE))
            }
            // SignalToNoiseEstimatorMedian.h:195-203. The source member is an
            // integer parameter in 0..=100, so it is never NaN.
            NoiseHistogramRange::Percentile { percentile }
                if !(0.0..=100.0).contains(&percentile) =>
            {
                Err(source_invalid_value(percentile, PERCENTILE_MESSAGE))
            }
            _ => Ok(()),
        }
    }

    /// Estimate signal-to-noise ratios for parallel position and intensity
    /// slices under the native safety profile.
    ///
    /// Source `init(container)` followed by `getSignalToNoise(i)` for every
    /// point, plus `getSparseWindowPercent` and `getHistogramRightmostPercent`.
    /// Positions must be finite, strictly increasing and in coordinate units;
    /// intensities finite and nonnegative. Use
    /// [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`] for the
    /// source's behaviour.
    ///
    /// The source's containers store `f32` intensities; values here that are
    /// not `f32` values have no source counterpart. The percentile range
    /// narrows each intensity to `f32` for its pre-histogram arithmetic, as the
    /// source computes it in `f32`.
    ///
    /// The standard-deviation range sums intensities and squared deviations in
    /// input order, as the source does. An empty input returns empty estimates,
    /// a zero upper end in the standard-deviation range and zero percentages;
    /// the source divides zero by zero there (see
    /// [`NoiseCompatibility::nan_for_empty_input`]).
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] for invalid options or resource limits, a
    ///   manual range that is not positive or a percentile outside `0..=100`
    ///   (the source's `Exception::InvalidValue`, with its text), mismatched
    ///   slice lengths, non-finite values, negative intensities or duplicate
    ///   positions, a non-finite histogram upper end, window bound or ratio, or
    ///   an exhausted work budget.
    /// * [`Error::UnsortedData`] when positions decrease.
    /// * [`Error::Unsupported`] where the source is undefined: more than
    ///   `i32::MAX` points in the manual and standard-deviation ranges, and in
    ///   the percentile range unless the range comes out negative (then only a
    ///   pre-histogram bin above `i32::MAX` points, `:216`, or a running count
    ///   above `i32::MAX`, `:228`); a window of exactly `i32::MAX` points
    ///   (`:324`); and the percentile range outside its domain (see
    ///   [`NoiseHistogramRange::Percentile`]).
    pub fn estimate(&self, positions: &[f64], intensities: &[f64]) -> Result<NoiseEstimates> {
        self.estimate_with_compatibility(positions, intensities, &PickingCompatibility::default())
    }

    /// As [`SignalToNoiseEstimatorMedian::estimate`], adopting the source
    /// behaviours `compatibility` selects.
    ///
    /// With `allow_duplicate_positions` and `allow_unsorted_positions`, equal
    /// and decreasing positions are accepted; the sliding window then moves its
    /// two ends exactly as the source's iterators do, never past the centre. With
    /// `allow_negative_intensities`, negative intensities are accepted; they
    /// fall into the first histogram bin, as the source's clamped bin index
    /// puts them. If the automatic range then comes out negative, the source
    /// logs a warning and returns before computing anything, leaving every
    /// signal-to-noise ratio at `0.0`; this returns exactly that, with `noise`
    /// set to `f64::INFINITY`, both percentages `0.0` and the warning in
    /// [`NoiseEstimates::log`]. [`PickingCompatibility::noise`] selects the
    /// estimator's own source behaviours ([`NoiseCompatibility`]).
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::estimate`], without the refusals that
    /// `compatibility` lifts.
    pub fn estimate_with_compatibility(
        &self,
        positions: &[f64],
        intensities: &[f64],
        compatibility: &PickingCompatibility,
    ) -> Result<NoiseEstimates> {
        if positions.len() != intensities.len() {
            return Err(bad("signal arrays differ in length or exceed point limit"));
        }
        self.estimate_points(
            positions.len(),
            |i| positions[i],
            |i| intensities[i],
            compatibility,
            None,
        )
    }

    /// As [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`],
    /// reporting progress to `progress` as the source's `ProgressLogger` base
    /// does.
    ///
    /// The source calls `startProgress(0, n, "noise estimation of data")`
    /// after the histogram range is known and the negative-range early return
    /// has not been taken (`:288`), `setProgress(k)` after the `k`-th window
    /// (`:367`) and `endProgress()` after the last (`:371`), for an empty input
    /// too. The other entry points report nothing, as a source object whose
    /// log type is `NONE`. If a native refusal stops the estimation after the
    /// report started, the report is ended before the error is returned, so the
    /// logger's nesting stays balanced; the source cannot fail there.
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`], and the
    /// errors of `progress`.
    pub fn estimate_with_progress(
        &self,
        positions: &[f64],
        intensities: &[f64],
        compatibility: &PickingCompatibility,
        progress: &mut ProgressLogger,
    ) -> Result<NoiseEstimates> {
        if positions.len() != intensities.len() {
            return Err(bad("signal arrays differ in length or exceed point limit"));
        }
        self.estimate_points(
            positions.len(),
            |i| positions[i],
            |i| intensities[i],
            compatibility,
            Some(progress),
        )
    }

    /// Source `init(spectrum)` for `SignalToNoiseEstimatorMedian<MSSpectrum>`:
    /// the estimates of a spectrum's peaks, with their `f32` intensities.
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`].
    pub fn estimate_spectrum(
        &self,
        spectrum: &MSSpectrum,
        compatibility: &PickingCompatibility,
    ) -> Result<NoiseEstimates> {
        self.estimate_peaks(&spectrum.peaks, compatibility, None)
    }

    /// Source `init(chromatogram)` for
    /// `SignalToNoiseEstimatorMedian<MSChromatogram>`, the instantiation
    /// `PeakPickerHiRes` and `PeakPickerChromatogram` use: positions are
    /// retention times in seconds.
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`].
    pub fn estimate_chromatogram(
        &self,
        chromatogram: &MSChromatogram,
        compatibility: &PickingCompatibility,
    ) -> Result<NoiseEstimates> {
        self.estimate_peaks(&chromatogram.peaks, compatibility, None)
    }

    /// Source `init(container)` over any container's points, optionally with a
    /// progress report (see
    /// [`SignalToNoiseEstimatorMedian::estimate_with_progress`]).
    ///
    /// # Errors
    ///
    /// As [`SignalToNoiseEstimatorMedian::estimate_with_progress`].
    pub fn estimate_peaks<P: NoisePoint>(
        &self,
        peaks: &[P],
        compatibility: &PickingCompatibility,
        progress: Option<&mut ProgressLogger>,
    ) -> Result<NoiseEstimates> {
        self.estimate_points(
            peaks.len(),
            |i| peaks[i].position(),
            |i| x86::widen(peaks[i].intensity()),
            compatibility,
            progress,
        )
    }

    /// Estimate over picker points without copying them into slices.
    pub(super) fn estimate_signal<P: SignalPoint>(
        &self,
        points: &[P],
        compatibility: &PickingCompatibility,
    ) -> Result<NoiseEstimates> {
        self.estimate_points(
            points.len(),
            |i| points[i].position(),
            |i| x86::widen(points[i].intensity()),
            compatibility,
            None,
        )
    }

    fn validate_points(
        &self,
        n: usize,
        position: &impl Fn(usize) -> f64,
        intensity: &impl Fn(usize) -> f64,
        compatibility: &PickingCompatibility,
    ) -> Result<()> {
        if n > self.max_points {
            return Err(bad("signal arrays differ in length or exceed point limit"));
        }
        let finite_only = !compatibility.noise.source_value_domain;
        for i in 0..n {
            let (x, y) = (position(i), intensity(i));
            if finite_only && (!x.is_finite() || !y.is_finite()) {
                return Err(bad(
                    "noise estimation requires finite coordinates and intensities; non-finite values need NoiseCompatibility::source_value_domain",
                ));
            }
            if y < 0.0 && !compatibility.allow_negative_intensities {
                return Err(bad(
                    "noise estimation requires nonnegative intensities; negative intensities need PickingCompatibility::allow_negative_intensities",
                ));
            }
        }
        if !compatibility.allow_unsorted_positions && (1..n).any(|i| position(i - 1) > position(i))
        {
            return Err(Error::UnsortedData);
        }
        if !compatibility.allow_duplicate_positions
            && (1..n).any(|i| position(i - 1) == position(i))
        {
            return Err(bad(
                "noise estimation requires distinct coordinates; duplicates need PickingCompatibility::allow_duplicate_positions",
            ));
        }
        Ok(())
    }

    /// The histogram upper end, the source member `max_intensity_` after
    /// `SignalToNoiseEstimatorMedian.h:185-245`.
    fn histogram_upper_end(
        &self,
        n: usize,
        intensity: &impl Fn(usize) -> f64,
        compatibility: &NoiseCompatibility,
    ) -> Result<f64> {
        match self.histogram_range {
            NoiseHistogramRange::Manual { max_intensity } => Ok(max_intensity),
            NoiseHistogramRange::StandardDeviation { factor } => {
                if n == 0 && !compatibility.nan_for_empty_input {
                    return Ok(0.0);
                }
                if n > SOURCE_MAX_POINTS {
                    return Err(Error::Unsupported(format!(
                        "SignalToNoiseEstimatorMedian is undefined for {n} points: `++size` at \
                         SignalToNoiseEstimator.h:123 overflows int beyond {SOURCE_MAX_POINTS}"
                    )));
                }
                // :188-189 with SignalToNoiseEstimator::estimate_; the Release
                // build computes `sqrt(v) * factor + mean`.
                let gauss = GaussianEstimate::of_indexed(n, intensity);
                Ok(x86::add(
                    x86::mul(x86::sqrt(gauss.variance), factor),
                    gauss.mean,
                ))
            }
            NoiseHistogramRange::Percentile { percentile } => {
                Self::percentile_upper_end(n, intensity, percentile)
            }
        }
    }

    /// The histogram upper end of AUTOMAXBYPERCENT
    /// (`SignalToNoiseEstimatorMedian.h:205-232`) on the source's defined
    /// domain; see [`NoiseHistogramRange::Percentile`].
    fn percentile_upper_end(
        n: usize,
        intensity: &impl Fn(usize) -> f64,
        percentile: f64,
    ) -> Result<f64> {
        if n == 0 {
            return Err(percentile_undefined(
                "SignalToNoiseEstimatorMedian.h:209 dereferences end() of an empty container",
            ));
        }
        // The source reads `float` intensities; `as f32` is the identity on them.
        let narrow = |i: usize| intensity(i) as f32;
        // :208-209, `minss`: the first minimum under `<`.
        let mut minimum = narrow(0);
        for i in 1..n {
            let value = narrow(i);
            if value < minimum {
                minimum = value;
            }
        }
        // :211, a float division widened to double.
        let bin_size = f64::from(minimum / 100.0_f32);
        let mut histogram = [0usize; 100];
        for i in 0..n {
            // :216, `subss 1.0f; cvtss2sd; divsd bin_size; cvttsd2si`.
            let quotient = f64::from(narrow(i) - 1.0_f32) / bin_size;
            if !(quotient > -1.0 && quotient < 100.0) {
                return Err(percentile_undefined(&format!(
                    "at point {i} (intensity {}), SignalToNoiseEstimatorMedian.h:216 writes the \
                     pre-histogram at the truncation of {quotient}, outside [0, 99] (minimum \
                     intensity {minimum}, bin size {bin_size})",
                    narrow(i)
                )));
            }
            // -1 < quotient < 100: the truncation is in 0..=99.
            count_in_pre_histogram(&mut histogram, quotient as usize, i)?;
        }
        let last = percentile_walk(n, &histogram, percentile)?;
        // :232, `(i + 0.5) * bin_size` with i = -1 when no bin was visited.
        let index = last.map_or(-1.0, |b| b as f64);
        Ok((index + 0.5) * bin_size)
    }

    fn estimate_points(
        &self,
        n: usize,
        position: impl Fn(usize) -> f64,
        intensity: impl Fn(usize) -> f64,
        compatibility: &PickingCompatibility,
        mut progress: Option<&mut ProgressLogger>,
    ) -> Result<NoiseEstimates> {
        let noise_compat = compatibility.noise;
        let source = noise_compat.source_value_domain;
        self.validate(&noise_compat)?;
        self.validate_points(n, &position, &intensity, compatibility)?;
        let max_intensity = self.histogram_upper_end(n, &intensity, &noise_compat)?;
        if !source && !max_intensity.is_finite() && n != 0 {
            return Err(bad(
                "noise histogram maximum is not finite; this needs NoiseCompatibility::source_value_domain",
            ));
        }
        let mut result = NoiseEstimates {
            signal_to_noise: Vec::new(),
            noise: Vec::new(),
            max_intensity,
            sparse_window_percent: 0.0,
            histogram_rightmost_percent: 0.0,
            log: Vec::new(),
        };
        if max_intensity < 0.0 {
            // :247-251: an unconditional warning, then a return with the
            // zero-initialised estimates.
            result.log.push(format!(
                "SignalToNoiseEstimatorMedian: the max_intensity_ value should be positive! {}",
                stream_double(max_intensity)
            ));
            result.signal_to_noise.resize(n, 0.0);
            result.noise.resize(n, f64::INFINITY);
            return Ok(result);
        }
        // Refused before the estimates are allocated.
        main_loop_counters_fit(n)?;
        result.signal_to_noise.reserve_exact(n);
        result.noise.reserve_exact(n);
        let mut histogram = vec![0usize; self.bin_count];
        if let Some(logger) = progress.as_deref_mut() {
            logger.start_progress(0, n as i64, NOISE_PROGRESS_LABEL)?;
        }
        let outcome = self.windows(
            n,
            &position,
            &intensity,
            max_intensity,
            &noise_compat,
            &mut histogram,
            &mut result,
            &mut progress,
        );
        let (sparse, rightmost) = match progress {
            Some(logger) => {
                let ended = logger.end_progress(0);
                let counts = outcome?;
                ended?;
                counts
            }
            None => outcome?,
        };
        if n == 0 {
            if noise_compat.nan_for_empty_input {
                let nan = f64::from_bits(x86::DEFAULT_NAN);
                result.sparse_window_percent = nan;
                result.histogram_rightmost_percent = nan;
            }
        } else {
            // :373-374, `(count * 100) / window_count`.
            result.sparse_window_percent = sparse * 100.0 / n as f64;
            result.histogram_rightmost_percent = rightmost * 100.0 / n as f64;
        }
        // :377-393; a NaN percentage compares false.
        if result.sparse_window_percent > 20.0 && self.write_log_messages {
            result.log.push(format!(
                "WARNING in SignalToNoiseEstimatorMedian: {}% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'",
                stream_double(result.sparse_window_percent)
            ));
        }
        if result.histogram_rightmost_percent > 1.0 && self.write_log_messages {
            result.log.push(format!(
                "WARNING in SignalToNoiseEstimatorMedian: {}% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)",
                stream_double(result.histogram_rightmost_percent)
            ));
        }
        Ok(result)
    }

    /// The main loop of `computeSTN_` (`:253-369`). Returns the sparse-window
    /// and rightmost-bin counts.
    #[allow(clippy::too_many_arguments)]
    fn windows(
        &self,
        n: usize,
        position: &impl Fn(usize) -> f64,
        intensity: &impl Fn(usize) -> f64,
        max_intensity: f64,
        compatibility: &NoiseCompatibility,
        histogram: &mut [usize],
        result: &mut NoiseEstimates,
        progress: &mut Option<&mut ProgressLogger>,
    ) -> Result<(f64, f64)> {
        let source = compatibility.source_value_domain;
        // :257-259, `0.5 * win_len` and `maxsd`, which maps a NaN to 1.
        let half_window = 0.5 * self.window_length;
        let width = x86::max_one(max_intensity / self.bin_count as f64);
        let last_bin = self.bin_count - 1;
        let conversion = compatibility.bin_index;
        let to_bin = |y: f64| bin_index(y, width, last_bin, conversion);
        let (mut left, mut right, mut work) = (0usize, 0usize, 0usize);
        let mut charge = || -> Result<()> {
            if work == self.max_work {
                Err(bad("noise estimation exceeds configured work limit"))
            } else {
                work += 1;
                Ok(())
            }
        };
        let mut sparse = 0.0;
        let mut rightmost = 0.0;
        for i in 0..n {
            let center = position(i);
            let low = center - half_window;
            let high = center + half_window;
            if !source && (!low.is_finite() || !high.is_finite()) {
                return Err(bad(
                    "noise window coordinates overflow; this needs NoiseCompatibility::source_value_domain",
                ));
            }
            // The left end never passes the centre, and it only removes points
            // the right end added, even for unsorted, duplicate or non-finite
            // positions (see the support document), so `left <= right` holds
            // and the guard never changes the source's iteration.
            while left < right && position(left) < low {
                charge()?;
                histogram[to_bin(intensity(left))] -= 1;
                left += 1;
            }
            while right < n && position(right) <= high {
                charge()?;
                histogram[to_bin(intensity(right))] += 1;
                right += 1;
            }
            let count = right - left;
            let noise = if count < self.min_required_elements {
                sparse += 1.0;
                self.noise_for_empty_window
            } else {
                if count == SOURCE_MAX_POINTS {
                    return Err(Error::Unsupported(
                        "SignalToNoiseEstimatorMedian is undefined here: \
                         `elements_in_window + 1` overflows int at SignalToNoiseEstimatorMedian.h:324"
                            .into(),
                    ));
                }
                // :322-329: the first bin whose cumulative count reaches
                // ceil(count / 2), never past the last bin. The source starts
                // from bin -1 with a zero count; the first iteration always runs
                // (count >= 1 and bin_count >= 3), which this unrolls.
                let half = count.div_ceil(2);
                let mut median = 0usize;
                let mut cumulative = histogram[0];
                while median < last_bin && cumulative < half {
                    charge()?;
                    median += 1;
                    cumulative += histogram[median];
                }
                if median == last_bin {
                    rightmost += 1.0;
                }
                let in_bin = histogram[median];
                // The bin the walk stops in is never empty: the walk entered it
                // with a cumulative count below `half` and every bin count is
                // the number of window points in that bin, which sum to `count`
                // >= `half`, so the stopping bin holds at least one point. The
                // source's fallback to the bin centre (`:350-353`, whose comment
                // says it runs when the rightmost bin is hit while empty) is
                // therefore unreachable; it is kept for fidelity.
                let estimate = if in_bin > 0 {
                    let before = cumulative - in_bin;
                    median as f64 * width + (half - before) as f64 / in_bin as f64 * width
                } else {
                    (median as f64 + 0.5) * width
                };
                // :356, `std::max(1.0, estimate)`.
                x86::max_one(estimate)
            };
            // :360, `cvtss2sd I; divsd noise`: the intensity is the first operand.
            let ratio = x86::div(intensity(i), noise);
            if !source && (!noise.is_finite() || !ratio.is_finite()) {
                return Err(bad(
                    "noise estimate is not finite; this needs NoiseCompatibility::source_value_domain",
                ));
            }
            result.noise.push(noise);
            result.signal_to_noise.push(ratio);
            if let Some(logger) = progress.as_deref_mut() {
                logger.set_progress(i as i64 + 1)?;
            }
        }
        Ok((sparse, rightmost))
    }
}

/// The refusal of an AUTOMAXBYPERCENT input on which the source is undefined.
fn percentile_undefined(detail: &str) -> Error {
    Error::Unsupported(format!(
        "SignalToNoiseEstimatorMedian auto_mode 1 (AUTOMAXBYPERCENT) is undefined for this input: {detail}"
    ))
}

/// `++histogram_auto[bin]` for the point `point`
/// (`SignalToNoiseEstimatorMedian.h:216`), whose counts are `int`: a bin that
/// already holds `i32::MAX` points overflows.
fn count_in_pre_histogram(histogram: &mut [usize; 100], bin: usize, point: usize) -> Result<()> {
    let slot = histogram.get_mut(bin).ok_or_else(|| {
        percentile_undefined("SignalToNoiseEstimatorMedian.h:216 writes outside the pre-histogram")
    })?;
    if *slot == SOURCE_MAX_POINTS {
        return Err(percentile_undefined(&format!(
            "at point {point}, `++histogram_auto[{bin}]` at SignalToNoiseEstimatorMedian.h:216 \
             overflows int beyond {SOURCE_MAX_POINTS}"
        )));
    }
    *slot += 1;
    Ok(())
}

/// The percentile walk of AUTOMAXBYPERCENT over `n` points
/// (`SignalToNoiseEstimatorMedian.h:220-230`): the last pre-histogram bin
/// visited, or `None` for the source's `i = -1`.
fn percentile_walk(n: usize, histogram: &[usize; 100], percentile: f64) -> Result<Option<usize>> {
    // :220, `(int)(auto_max_percentile_ * c.size() / 100)`: `cvtsi2sd n;
    // mulsd p; divsd 100.0; cvttsd2si` into a 32-bit register, then `test;
    // jle`, which skips the walk for a target <= 0. The product is never
    // negative or NaN (0 <= p <= 100 was checked); at or above 2^31 it does not
    // fit int, which is undefined behaviour, and the Release build's 32-bit
    // `cvttsd2si` returns INT_MIN there, which skips the walk. That is the
    // measured Linux x86-64 Release outcome (`sne-fix` oracle: n = 2^31 with
    // p = 100, and n = 3,000,000,001 with p = 72), reproduced here.
    let target = x86::cvttsd2si32(n as f64 * percentile / 100.0);
    let target = usize::try_from(target).unwrap_or(0);
    // :221-230: `elements_seen < elements_below_percentile` on int, at most one
    // bin per container element. A target <= 0 is 0 here: `seen` is never
    // negative, so the comparison is the same.
    let mut seen = 0usize;
    let mut last: Option<usize> = None;
    let mut run = 0usize;
    while run != n && seen < target {
        let bin = last.map_or(0, |b| b + 1);
        // target <= n (rounding is monotone, and 100 * n / 100 is exact for
        // any n a container can hold), and the 100 bins hold all n points, so
        // the walk stops by bin 99; `get` keeps the proof checked.
        let count = *histogram.get(bin).ok_or_else(|| {
            percentile_undefined(
                "the pre-histogram walk at SignalToNoiseEstimatorMedian.h:228 passes bin 99",
            )
        })?;
        // :228, `elements_seen += histogram_auto[i]` on int.
        if count > SOURCE_MAX_POINTS - seen {
            return Err(percentile_undefined(&format!(
                "`elements_seen += histogram_auto[{bin}]` at SignalToNoiseEstimatorMedian.h:228 \
                 overflows int ({seen} + {count})"
            )));
        }
        seen += count;
        last = Some(bin);
        run += 1;
    }
    Ok(last)
}

/// The main loop of `computeSTN_` counts windows (`++window_count`,
/// `SignalToNoiseEstimatorMedian.h:365`) and window elements
/// (`++elements_in_window`, `:310`) in `int`: beyond `i32::MAX` points the
/// window count overflows at the latest.
fn main_loop_counters_fit(n: usize) -> Result<()> {
    if n > SOURCE_MAX_POINTS {
        return Err(Error::Unsupported(format!(
            "SignalToNoiseEstimatorMedian is undefined for {n} points: its int counters overflow \
             beyond {SOURCE_MAX_POINTS}, `++window_count` at SignalToNoiseEstimatorMedian.h:365 at \
             the latest and `++elements_in_window` at :310 first if a window holds more points"
        )));
    }
    Ok(())
}

/// The histogram bin of an intensity (`SignalToNoiseEstimatorMedian.h:297`,
/// `:308`); see [`BinIndexConversion`].
fn bin_index(intensity: f64, width: f64, last_bin: usize, conversion: BinIndexConversion) -> usize {
    let quotient = intensity / width;
    match conversion {
        // `as usize` saturates a negative quotient to zero; `f64::min` returns
        // the bound for a NaN quotient.
        BinIndexConversion::ClampBeforeTruncation => quotient.min(last_bin as f64) as usize,
        BinIndexConversion::X86_64Release => {
            // bin_count <= i32::MAX was checked, so the last bin fits i32.
            let last = i32::try_from(last_bin).unwrap_or(i32::MAX);
            let truncated = x86::cvttsd2si32(quotient);
            usize::try_from(truncated.min(last).max(0)).unwrap_or(0)
        }
    }
}

impl SignalToNoiseEstimator for SignalToNoiseEstimatorMedian {
    type Estimates = NoiseEstimates;

    /// [`SignalToNoiseEstimatorMedian::estimate`], the native safety profile;
    /// [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`] selects
    /// the source behaviour.
    fn compute_stn(&self, positions: &[f64], intensities: &[f64]) -> Result<NoiseEstimates> {
        self.estimate(positions, intensities)
    }

    fn signal_to_noise(estimates: &NoiseEstimates) -> &[f64] {
        &estimates.signal_to_noise
    }

    /// [`SignalToNoiseEstimatorMedian::estimate_spectrum`] with the native
    /// safety profile.
    fn compute_stn_spectrum(&self, spectrum: &MSSpectrum) -> Result<NoiseEstimates> {
        self.estimate_spectrum(spectrum, &PickingCompatibility::default())
    }

    /// [`SignalToNoiseEstimatorMedian::estimate_chromatogram`] with the native
    /// safety profile.
    fn compute_stn_chromatogram(&self, chromatogram: &MSChromatogram) -> Result<NoiseEstimates> {
        self.estimate_chromatogram(chromatogram, &PickingCompatibility::default())
    }
}

#[cfg(test)]
mod tests {
    //! The `int` counter refusals beyond `i32::MAX` points, on the helpers that
    //! take the point count, since a full estimation there needs tens of
    //! gigabytes. `../oracle/sne-fix` runs the full estimation at `n = 2^31`
    //! against the Linux x86-64 Release build.
    use super::*;
    use std::cell::Cell;

    const TWO_30: usize = 1 << 30;
    const TWO_31: usize = 1 << 31;

    /// The pre-histogram of the Release cases: intensities alternating `2.0`
    /// and `2.5` fall into bins 50 and 75 (quotients `50.0000011...` and
    /// `75.0000016...` with `bin_size = 2.0f / 100.0f`).
    fn two_bins(first: usize, second: usize) -> [usize; 100] {
        let mut histogram = [0; 100];
        histogram[50] = first;
        histogram[75] = second;
        histogram
    }

    fn unsupported_at(result: Result<impl std::fmt::Debug>, line: &str) {
        match result {
            Err(Error::Unsupported(message)) => {
                assert!(message.contains(line), "{message} does not name {line}");
            }
            other => panic!("expected Unsupported naming {line}, got {other:?}"),
        }
    }

    #[test]
    fn percentile_walk_is_defined_beyond_int_max_points_when_the_target_skips_it() -> Result<()> {
        let histogram = two_bins(TWO_30, TWO_30);
        // :220, 0 * 2^31 / 100 = 0: no walk, i = -1 (Release case
        // chromatogram_2147483648_p0 and spectrum_2147483648_p0).
        assert_eq!(percentile_walk(TWO_31, &histogram, 0.0)?, None);
        // 100 * 2^31 / 100 = 2^31 does not fit int: the 32-bit cvttsd2si gives
        // INT_MIN and `test; jle` skips the walk (Release case *_p100).
        assert_eq!(percentile_walk(TWO_31, &histogram, 100.0)?, None);
        // 72 * 3,000,000,001 / 100 = 2,160,000,000.72 >= 2^31 (Release case
        // chromatogram_3000000001_p72): INT_MIN again.
        let big = two_bins(1_500_000_001, 1_500_000_000);
        assert_eq!(percentile_walk(3_000_000_001, &big, 72.0)?, None);
        // A tiny percentile truncates to zero: 4e-8 * 2^31 / 100 = 0.859.
        assert_eq!(percentile_walk(TWO_31, &histogram, 4e-8)?, None);
        Ok(())
    }

    #[test]
    fn percentile_walk_refuses_exactly_where_its_running_count_overflows() -> Result<()> {
        // 49 * 2^31 / 100 = 1,052,266,987.52 <= 2^30: the walk stops in bin 50,
        // with a positive upper end, so the main loop refuses the point count.
        let histogram = two_bins(TWO_30, TWO_30);
        assert_eq!(percentile_walk(TWO_31, &histogram, 49.0)?, Some(50));
        unsupported_at(
            main_loop_counters_fit(TWO_31),
            "SignalToNoiseEstimatorMedian.h:365",
        );
        // 99 * 2^31 / 100 = 2,126,008,811.52 > 2^30: bin 75 takes the count
        // to 2^31, one past INT_MAX (:228).
        unsupported_at(
            percentile_walk(TWO_31, &histogram, 99.0),
            "SignalToNoiseEstimatorMedian.h:228",
        );
        // One point fewer ends the walk exactly at INT_MAX, which int holds:
        // 100 * (2^31 - 1) / 100 = 2^31 - 1 fits the conversion too.
        let at_max = two_bins(TWO_30, TWO_30 - 1);
        assert_eq!(
            percentile_walk(SOURCE_MAX_POINTS, &at_max, 100.0)?,
            Some(75)
        );
        main_loop_counters_fit(SOURCE_MAX_POINTS)?;
        Ok(())
    }

    #[test]
    fn pre_histogram_counts_refuse_exactly_at_int_max() -> Result<()> {
        let mut histogram = [0; 100];
        histogram[7] = SOURCE_MAX_POINTS - 1;
        count_in_pre_histogram(&mut histogram, 7, 0)?;
        assert_eq!(histogram[7], SOURCE_MAX_POINTS);
        unsupported_at(
            count_in_pre_histogram(&mut histogram, 7, 1),
            "SignalToNoiseEstimatorMedian.h:216",
        );
        assert_eq!(histogram[7], SOURCE_MAX_POINTS);
        // Other bins are independent.
        count_in_pre_histogram(&mut histogram, 8, 2)?;
        Ok(())
    }

    #[test]
    fn standard_deviation_range_refuses_beyond_int_max_points_before_reading() {
        let reads = Cell::new(0usize);
        let intensity = |_: usize| {
            reads.set(reads.get() + 1);
            1.0
        };
        let estimator = SignalToNoiseEstimatorMedian::default();
        unsupported_at(
            estimator.histogram_upper_end(TWO_31, &intensity, &NoiseCompatibility::source()),
            "SignalToNoiseEstimator.h:123",
        );
        assert_eq!(reads.get(), 0);
        // The manual range has no counter of its own there.
        let manual = SignalToNoiseEstimatorMedian {
            histogram_range: NoiseHistogramRange::Manual { max_intensity: 5.0 },
            ..SignalToNoiseEstimatorMedian::default()
        };
        assert_eq!(
            manual
                .histogram_upper_end(TWO_31, &intensity, &NoiseCompatibility::source())
                .ok(),
            Some(5.0)
        );
        assert_eq!(reads.get(), 0);
    }
}
