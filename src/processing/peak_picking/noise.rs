// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Histogram-median signal-to-noise estimation.
//!
//! Ports the header-only template
//! `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h`
//! and the parts of its base `SignalToNoiseEstimator.h` it uses (the Gaussian
//! mean and variance estimate). `docs/PEAK_PICKING_SUPPORT.md` holds the API
//! mapping, the preserved conventions, the native differences and the evidence.
//!
//! The C++ class is a `DefaultParamHandler` that computes its estimates in
//! `init` and serves them through `getSignalToNoise(index)`. This port returns
//! all estimates at once from
//! [`SignalToNoiseEstimatorMedian::estimate`](crate::processing::peak_picking::SignalToNoiseEstimatorMedian::estimate),
//! and maps the parameter contract through `defaults`, `from_param` and
//! `to_param` on the same type.

use super::{PickingCompatibility, SignalPoint, bad};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};

/// How the upper end of the intensity histogram is chosen.
///
/// Source `SignalToNoiseEstimatorMedian::IntensityThresholdCalculation` with the
/// parameter `auto_mode`: `MANUAL = -1`, `AUTOMAXBYSTDEV = 0` (the default) and
/// `AUTOMAXBYPERCENT = 1`. Intensities at or above the upper end go into the
/// last histogram bin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NoiseHistogramRange {
    /// `auto_mode = 0`: mean plus `factor` times the population standard
    /// deviation of all intensities (`auto_max_stdev_factor`, source range
    /// `0..=999`).
    StandardDeviation {
        /// Multiplier of the standard deviation.
        factor: f64,
    },
    /// `auto_mode = -1`: a fixed upper end (`max_intensity`). The source throws
    /// `Exception::InvalidValue` when estimation runs with a value `<= 0`; this
    /// returns [`Error::InvalidValue`] at the same point.
    Manual {
        /// Upper end of the histogram. The source parameter is an integer; any
        /// positive finite value is accepted here.
        max_intensity: f64,
    },
    /// `auto_mode = 1`: the `auto_max_percentile`-th percentile of a
    /// 100-bin pre-histogram.
    ///
    /// Not implemented: at the pinned revision the source finds the *minimum*
    /// intensity with a reversed `std::max_element` comparator, divides by a
    /// bin size derived from it and indexes its pre-histogram without a bounds
    /// check. The executed product SDK crashes on ordinary data (`SIGBUS` or
    /// `SIGSEGV`). [`SignalToNoiseEstimatorMedian::estimate`] therefore returns
    /// [`Error::Unsupported`] whenever estimation runs in this mode. Selecting
    /// the mode is accepted, because the source only crashes when the estimator
    /// is initialised: a picker with `signal_to_noise = 0` never is.
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

/// Estimates the signal-to-noise ratio of every data point from the median of
/// a sliding-window intensity histogram.
///
/// For each data point, the points within half of
/// [`window_length`](Self::window_length) on either side (both ends inclusive)
/// form its window. The noise is the median intensity of that window, read
/// from a histogram of [`bin_count`](Self::bin_count) bins over
/// `[0, upper end]`: the median bin is the first whose cumulative count reaches
/// `(count + 1) / 2`, and the noise is linearly interpolated inside that bin,
/// assuming the bin's points are uniformly distributed, then floored at one.
/// Bins are at least one intensity unit wide. A window with fewer than
/// [`min_required_elements`](Self::min_required_elements) points is *sparse*
/// and uses [`noise_for_empty_window`](Self::noise_for_empty_window). The
/// signal-to-noise ratio is the point's intensity divided by its noise.
///
/// # Notes
///
/// The source warns through `OPENMS_LOG_WARN` when more than 20 percent of
/// windows are sparse or more than 1 percent of medians fall into the last bin
/// (an unreliable estimate; increase `max_intensity` and possibly `bin_count`),
/// unless `write_log_messages` is off. This port does not log; both percentages
/// are returned in [`NoiseEstimates`], and
/// [`write_log_messages`](Self::write_log_messages) is kept only as a parameter
/// value.
///
/// The source recomputes lazily when a parameter changes; this port recomputes
/// on every call and has no cache to invalidate.
#[derive(Clone, Debug, PartialEq)]
pub struct SignalToNoiseEstimatorMedian {
    /// Upper end of the intensity histogram (`auto_mode` and the parameter it
    /// selects).
    pub histogram_range: NoiseHistogramRange,
    /// Values of the range parameters the selected mode does not use; see
    /// [`NoiseRangeParameters`].
    pub range_parameters: NoiseRangeParameters,
    /// Full window width `win_len` in coordinate units (Th for spectra,
    /// seconds for chromatograms). At least one; source default `200`.
    pub window_length: f64,
    /// Number of histogram bins `bin_count`, at least three; source default
    /// `30`.
    pub bin_count: usize,
    /// Minimum points in a window, `min_required_elements`, at least one;
    /// source default `10`.
    pub min_required_elements: usize,
    /// Noise used for sparse windows, `noise_for_empty_window`; source default
    /// `1e20`. Must be finite and positive here; the source accepts any value
    /// (see the support document).
    pub noise_for_empty_window: f64,
    /// Source `write_log_messages`. Retained for the parameter contract; this
    /// port never logs.
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
    /// Signal-to-noise ratio per input point, source `getSignalToNoise(i)`.
    pub signal_to_noise: Vec<f64>,
    /// Noise per input point (native). `f64::INFINITY` for every point when the
    /// source returns early; see [`SignalToNoiseEstimatorMedian::estimate`].
    pub noise: Vec<f64>,
    /// The histogram upper end that was used.
    pub max_intensity: f64,
    /// Percentage of sparse windows, source `getSparseWindowPercent`.
    pub sparse_window_percent: f64,
    /// Percentage of medians found in the last bin, source
    /// `getHistogramRightmostPercent`.
    pub histogram_rightmost_percent: f64,
}

/// Parameter handler name of the source class, used in parameter diagnostics.
pub const SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME: &str = "SignalToNoiseEstimatorMedian";

const MAX_INTENSITY_DESCRIPTION: &str = "maximal intensity considered for histogram construction. By default, it will be calculated automatically (see auto_mode). Only provide this parameter if you know what you are doing (and change 'auto_mode' to '-1')! All intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram bin. If you choose 'max_intensity' too small, the noise estimate might be too small as well.  If chosen too big, the bins become quite large (which you could counter by increasing 'bin_count', which increases runtime). In general, the Median-S/N estimator is more robust to a manual max_intensity than the MeanIterative-S/N.";

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

impl SignalToNoiseEstimatorMedian {
    /// The source parameter defaults, as `SignalToNoiseEstimatorMedian()` sets
    /// them: names, values, value types, descriptions, restrictions and the
    /// `advanced` tag, in declaration order.
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
    /// Source `getParameters`.
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
        let mut handler = DefaultParamHandler::new(SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
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

    pub(super) fn validate(&self) -> Result<()> {
        if !self.window_length.is_finite()
            || self.window_length < 1.0
            || self.bin_count < 3
            || self.bin_count > self.max_bins
            || self.min_required_elements == 0
            || !self.noise_for_empty_window.is_finite()
            || self.noise_for_empty_window <= 0.0
            || self.max_points == 0
            || self.max_work == 0
        {
            return Err(bad("invalid median-noise options or resource limits"));
        }
        match self.histogram_range {
            NoiseHistogramRange::StandardDeviation { factor }
                if !factor.is_finite() || !(0.0..=999.0).contains(&factor) =>
            {
                Err(bad("noise standard-deviation factor must be in 0..=999"))
            }
            NoiseHistogramRange::Manual { max_intensity }
                if !max_intensity.is_finite() || max_intensity <= 0.0 =>
            {
                Err(bad(
                    "auto_mode is on MANUAL! max_intensity is <=0. Needs to be positive!",
                ))
            }
            NoiseHistogramRange::Percentile { percentile }
                if !percentile.is_finite() || !(0.0..=100.0).contains(&percentile) =>
            {
                Err(bad(
                    "auto_mode is on AUTOMAXBYPERCENT! auto_max_percentile is not in [0,100].",
                ))
            }
            NoiseHistogramRange::Percentile { .. } => Err(Error::Unsupported(
                "SignalToNoiseEstimatorMedian auto_mode 1 (AUTOMAXBYPERCENT) reads out of bounds in the source and is not ported".into(),
            )),
            _ => Ok(()),
        }
    }

    /// Estimate signal-to-noise ratios for parallel position and intensity
    /// slices.
    ///
    /// Source `init(container)` followed by `getSignalToNoise(i)` for every
    /// point, plus `getSparseWindowPercent` and `getHistogramRightmostPercent`.
    /// Positions must be finite, strictly increasing and in coordinate units;
    /// intensities finite and nonnegative. Use
    /// [`SignalToNoiseEstimatorMedian::estimate_with_compatibility`] for the
    /// source's acceptance of duplicate positions and negative intensities.
    ///
    /// The standard-deviation range sums intensities and squared deviations in
    /// input order, as the source does. An empty input returns empty estimates
    /// and zero percentages (the source divides zero by zero).
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidValue`] for invalid options or resource limits, a
    ///   manual range that is not positive (the source's
    ///   `Exception::InvalidValue`), mismatched slice lengths, non-finite
    ///   values, negative intensities or duplicate positions, or an exhausted
    ///   work budget.
    /// * [`Error::UnsortedData`] when positions decrease.
    /// * [`Error::Unsupported`] for [`NoiseHistogramRange::Percentile`].
    pub fn estimate(&self, positions: &[f64], intensities: &[f64]) -> Result<NoiseEstimates> {
        self.estimate_with_compatibility(positions, intensities, &PickingCompatibility::default())
    }

    /// As [`SignalToNoiseEstimatorMedian::estimate`], accepting the inputs the
    /// source accepts when `compatibility` says so.
    ///
    /// With `allow_duplicate_positions` and `allow_unsorted_positions`, equal
    /// and decreasing positions are accepted; the sliding window then moves its
    /// two ends exactly as the source's iterators do, never past the centre. With
    /// `allow_negative_intensities`, negative intensities are accepted; they
    /// fall into the first histogram bin, as the source's clamped bin index
    /// puts them. If the standard-deviation range then comes out negative, the
    /// source logs a warning and returns before computing anything, leaving
    /// every signal-to-noise ratio at `0.0`; this returns exactly that, with
    /// `noise` set to `f64::INFINITY` and both percentages `0.0`.
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
            |i| f64::from(points[i].intensity()),
            compatibility,
        )
    }

    fn estimate_points(
        &self,
        n: usize,
        position: impl Fn(usize) -> f64,
        intensity: impl Fn(usize) -> f64,
        compatibility: &PickingCompatibility,
    ) -> Result<NoiseEstimates> {
        self.validate()?;
        super::validate_points(n, &position, &intensity, self.max_points, compatibility)?;
        let max_intensity = match self.histogram_range {
            NoiseHistogramRange::Manual { max_intensity } => max_intensity,
            NoiseHistogramRange::StandardDeviation { factor } => {
                if n == 0 {
                    0.0
                } else {
                    // SignalToNoiseEstimator::estimate_: two in-order passes,
                    // population variance.
                    let mut mean = 0.0;
                    for i in 0..n {
                        mean += intensity(i);
                    }
                    mean /= n as f64;
                    let mut variance = 0.0;
                    for i in 0..n {
                        let deviation = mean - intensity(i);
                        variance += deviation * deviation;
                    }
                    variance /= n as f64;
                    mean + variance.sqrt() * factor
                }
            }
            // validate() refuses this mode before estimation.
            NoiseHistogramRange::Percentile { .. } => {
                return Err(Error::Unsupported(
                    "SignalToNoiseEstimatorMedian auto_mode 1 is not ported".into(),
                ));
            }
        };
        if !max_intensity.is_finite() {
            return Err(bad("noise histogram maximum overflows"));
        }
        let mut result = NoiseEstimates {
            signal_to_noise: Vec::with_capacity(n),
            noise: Vec::with_capacity(n),
            max_intensity,
            sparse_window_percent: 0.0,
            histogram_rightmost_percent: 0.0,
        };
        if max_intensity < 0.0 {
            // Source: warn "the max_intensity_ value should be positive!" and
            // return with the zero-initialised estimates.
            result.signal_to_noise.resize(n, 0.0);
            result.noise.resize(n, f64::INFINITY);
            return Ok(result);
        }
        let half_window = self.window_length / 2.0;
        let width = (max_intensity / self.bin_count as f64).max(1.0);
        let last_bin = self.bin_count - 1;
        // Source: clamp the truncated quotient to [0, bin_count - 1]; `as usize`
        // saturates negative quotients to zero.
        let to_bin = |y: f64| (y / width).min(last_bin as f64) as usize;
        let mut histogram = vec![0usize; self.bin_count];
        let (mut left, mut right, mut work) = (0, 0, 0usize);
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
            if !low.is_finite() || !high.is_finite() {
                return Err(bad("noise window coordinates overflow"));
            }
            // The centre never leaves its own window, so `left <= i < right`
            // holds after these loops; the `left < right` guard only bounds them.
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
                // Source loop: the first bin whose cumulative count reaches
                // ceil(count / 2), never past the last bin.
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
                let estimate = if in_bin > 0 {
                    let before = cumulative - in_bin;
                    median as f64 * width + (half - before) as f64 / in_bin as f64 * width
                } else {
                    // Only reachable when the last bin is hit while empty.
                    (median as f64 + 0.5) * width
                };
                estimate.max(1.0)
            };
            let ratio = intensity(i) / noise;
            if !noise.is_finite() || !ratio.is_finite() {
                return Err(bad("noise estimate is not finite"));
            }
            result.noise.push(noise);
            result.signal_to_noise.push(ratio);
        }
        if n != 0 {
            result.sparse_window_percent = sparse * 100.0 / n as f64;
            result.histogram_rightmost_percent = rightmost * 100.0 / n as f64;
        }
        Ok(result)
    }
}
