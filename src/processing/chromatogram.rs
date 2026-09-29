// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! OpenSWATH chromatogram picking with smoothed spline seeds, source boundary
//! extension and inclusive raw intensity sums. See CHROMATOGRAM_PICKING_SUPPORT.md.
//!
//! [`PeakPickerChromatogram::compatibility`] selects the source behaviours the
//! native default refuses; the internal noise estimate always reproduces the
//! source, because the signal it reads is this picker's own smoothed trace or a
//! chromatogram the picker has already validated.
//!
//! [`PeakPickerChromatogram::compatibility`]: crate::processing::chromatogram::PeakPickerChromatogram::compatibility

use super::checked_intensity;
use super::peak_picking::{
    FwhmUnit, PeakPickerHiRes, PickedChromatogram, PickingCompatibility,
    SignalToNoiseEstimatorMedian, flag, flag_entry, float, float_entry, int_entry, integer,
    string_entry,
};
use super::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use crate::kernel::{ChromatogramPeak, DataArray, MSChromatogram};
use crate::param::{DefaultParamHandler, Param, ParamNode, ParamValue};
use crate::{Error, Result};

/// The source class name, `DefaultParamHandler("PeakPickerChromatogram")`,
/// under which unknown-parameter warnings are reported.
pub const PEAK_PICKER_CHROMATOGRAM_NAME: &str = "PeakPickerChromatogram";

/// The source's exception text when `method` is `crawdad` in a build without
/// `WITH_CRAWDAD` (`PeakPickerChromatogram.cpp:414-420`), the default build.
const CRAWDAD_UNAVAILABLE: &str =
    "PeakPickerChromatogram was not compiled with crawdad, please choose a different algorithm!";

/// The source `method` parameter: which signal the peak boundaries are found on.
///
/// The source's third value, `crawdad`, wraps an external library that the
/// default build leaves out (`WITH_CRAWDAD` is off); it has no variant here, and
/// [`PeakPickerChromatogram::from_param`] refuses it with the source's own
/// exception text, as that build's `updateMembers_` does.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ChromatogramPickingMethod {
    /// Find boundaries on the original, unsmoothed signal.
    Legacy,
    /// Find boundaries on the smoothed signal; integrate original intensities.
    #[default]
    Corrected,
}

impl ChromatogramPickingMethod {
    /// The source parameter value: `legacy` or `corrected`.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Legacy => "legacy",
            Self::Corrected => "corrected",
        }
    }
}

/// The smoother applied before seeding: the source `use_gauss` switch with the
/// parameters of the selected filter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ChromatogramSmoothing {
    /// Width in seconds, equal to eight Gaussian standard deviations.
    Gaussian {
        /// Source `gauss_width`: Gaussian width in seconds, the estimated peak
        /// size.
        width: f64,
    },
    /// Fits sample indices; even frame lengths are incremented, as in OpenMS.
    SavitzkyGolay {
        /// Source `sgolay_frame_length`: the number of subsequent data points
        /// used for smoothing. This number has to be uneven; if it is not, one
        /// is added.
        frame_length: usize,
        /// Source `sgolay_polynomial_order`: order of the polynomial that is
        /// fitted.
        polynomial_order: usize,
    },
}
impl Default for ChromatogramSmoothing {
    fn default() -> Self {
        Self::Gaussian { width: 50.0 }
    }
}

/// Indices into the original input. Both boundary indices are inclusive.
/// Access original peaks to recover their full f64 RTs; the source metadata
/// arrays leftWidth/rightWidth store rounded f32 RT values.
///
/// Native: the source keeps these as the protected `left_width_` and
/// `right_width_` members and exposes only their `f32` retention times.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChromatogramPeakRegion {
    /// Sample closest to the seed apex (source `findClosestPeak_`).
    pub apex_index: usize,
    /// First sample of the peak (source `left_width_[i]`).
    pub left_index: usize,
    /// Last sample of the peak, inclusive (source `right_width_[i]`).
    pub right_index: usize,
}

/// What one [`PeakPickerChromatogram::pick_chromatogram`] call produces: the
/// source's `picked_chrom` and `smoothed_chrom` out-parameters and the native
/// integration regions.
#[derive(Clone, Debug, PartialEq)]
pub struct ChromatogramPickingResult {
    /// Spline seeds and five source arrays. `boundaries` are the HiRes seed
    /// support, while `regions` below holds the final integration boundaries.
    /// `omitted_arrays` explicitly lists unaggregated input profile arrays.
    pub picked: PickedChromatogram,
    /// Full original sampling and annotations, with smoothed intensities.
    pub smoothed: MSChromatogram,
    /// One integration region per picked peak, in picked order.
    pub regions: Vec<ChromatogramPeakRegion>,
    /// The boundary noise estimator's warnings (sparse windows, a median in
    /// the rightmost histogram bin), returned only when
    /// [`SignalToNoiseEstimatorMedian::write_log_messages`] is set, which the
    /// source's `write_sn_log_messages` parameter controls. The source writes
    /// them to its log stream; they are returned here. Empty when no noise
    /// estimate ran.
    pub noise_log: Vec<String>,
}

/// Finds peaks in a single chromatogram (source `PeakPickerChromatogram`).
///
/// It uses the [`PeakPickerHiRes`] internally to find interesting seed
/// candidates. These candidates are then expanded and a right/left border of
/// the peak is searched. Additionally, overlapping peaks can be removed.
///
/// The picked chromatogram's float data arrays report, per peak (indexed as
/// [`IDX_FWHM`](Self::IDX_FWHM) to [`IDX_SN`](Self::IDX_SN), source enum
/// `FLOATINDICES`): `FWHM`, `IntegratedIntensity`, `leftWidth`, `rightWidth`
/// and `SN` (the apex signal-to-noise ratio, as estimated by the internal
/// [`SignalToNoiseEstimatorMedian`]). `SN` is only populated with real values
/// if [`report_sn`](Self::report_sn) is enabled or
/// [`signal_to_noise`](Self::signal_to_noise) is above zero; otherwise it is
/// filled with -1 for every peak.
///
/// Native legacy/corrected methods only. Crawdad requires an external backend
/// and is not represented by either native method.
///
/// The source class is a `DefaultParamHandler`; its eleven parameters are
/// [`PeakPickerChromatogram::defaults`], [`PeakPickerChromatogram::to_param`]
/// and [`PeakPickerChromatogram::from_param`] (source `getDefaults`,
/// `getParameters`, and `setParameters` with `updateMembers_`). The source's
/// protected copy assignment operator is declared in the header but never
/// defined; [`Clone`] copies every field.
#[derive(Clone, Debug)]
pub struct PeakPickerChromatogram {
    /// Source `method`: `legacy` or `corrected` (default).
    pub method: ChromatogramPickingMethod,
    /// Source `use_gauss` (default `true`) with `gauss_width` (default 50 s)
    /// or `sgolay_frame_length` (default 15) and `sgolay_polynomial_order`
    /// (default 3).
    pub smoothing: ChromatogramSmoothing,
    /// The parameters of the smoother that [`smoothing`](Self::smoothing) does
    /// not select, kept so that [`to_param`](Self::to_param) writes back the
    /// values [`from_param`](Self::from_param) read: the source keeps all four
    /// smoothing parameters whichever filter `use_gauss` selects. Picking never
    /// reads this field. Default: Savitzky-Golay with frame 15 and order 3.
    pub inactive_smoothing: ChromatogramSmoothing,
    /// Extend through nondecreasing intensities inside this RT distance from
    /// the spline apex. The S/N check still applies. None disables this rule.
    pub peak_width: Option<f64>,
    /// Boundary-extension S/N threshold. Zero disables the boundary S/N gate.
    pub signal_to_noise: f64,
    /// Boundary and apex-report estimator (source `snt_`); the source window
    /// defaults to 1000 s (`sn_win_len`), the bin count to 30
    /// (`sn_bin_count`) and its log messages to off
    /// (`write_sn_log_messages`).
    ///
    /// This estimator always runs under [`PickingCompatibility::source`], so
    /// its whole source profile applies whatever
    /// [`compatibility`](Self::compatibility) says, including the Linux x86-64
    /// Release build's bin-index conversion: a histogram quotient that leaves
    /// `int` range is binned as that build bins it (bin `0`, from
    /// `cvttsd2si`'s `INT_MIN`) rather than clamped to the last bin. Producing
    /// one needs [`histogram_range`](SignalToNoiseEstimatorMedian::histogram_range)
    /// to be set by hand: the source's own picker sets only `win_len` and
    /// `bin_count` on `snt_`
    /// (`ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:408-412`), leaving
    /// `max_intensity` at `-1`, and with the automatic range the bin width
    /// scales with the data, so the quotient stays inside `int` for any
    /// histogram that fits in memory. This field exposes the whole estimator,
    /// which the source's parameter set does not.
    pub noise_estimator: SignalToNoiseEstimatorMedian,
    /// Seed picking is independent: C++ configures HiRes once at S/N 1, with
    /// its 200 s default noise window, even after boundary settings change.
    /// The source's constructor copies `signal_to_noise` into the seed picker
    /// before any `setParameters` call and `updateMembers_` never updates it,
    /// so [`from_param`](Self::from_param) leaves this at one.
    pub seed_signal_to_noise: f64,
    /// Noise estimator of the seed [`PeakPickerHiRes`] (source `pp_`), at that
    /// picker's defaults; no source parameter reaches it.
    pub seed_noise_estimator: SignalToNoiseEstimatorMedian,
    /// Source `report_sn`: report apex S/N even when the boundary S/N gate is
    /// disabled. Forces the noise estimate to run.
    pub report_sn: bool,
    /// Source `remove_overlapping_peaks`: try to remove overlapping peaks
    /// during peak picking.
    pub remove_overlapping_peaks: bool,
    /// Which source behaviours this picker adopts where the native default
    /// refuses; see [`PickingCompatibility`].
    ///
    /// Two flags change what this picker accepts.
    /// [`allow_duplicate_positions`](PickingCompatibility::allow_duplicate_positions)
    /// accepts equal retention times, which `MSChromatogram::isSorted` lets
    /// through to `snt_.init`, and
    /// [`allow_negative_intensities`](PickingCompatibility::allow_negative_intensities)
    /// accepts negative sample intensities. The flag is also handed to the
    /// seed [`PeakPickerHiRes`], as the source's `pp_` member reproduces the
    /// source unconditionally.
    ///
    /// [`allow_unsorted_positions`](PickingCompatibility::allow_unsorted_positions)
    /// has no effect here: `pickChromatogram`
    /// (`ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:68-72`) throws
    /// `Exception::IllegalArgument` for a chromatogram that is not sorted by
    /// position, so decreasing retention times stay refused in both profiles.
    ///
    /// The internal noise estimate does not consult this field; see
    /// [`PeakPickerChromatogram::pick_chromatogram`].
    pub compatibility: PickingCompatibility,
    /// Native ceiling on the input point count, checked before any work.
    pub max_points: usize,
    /// Per-stage work limit: conservative smoothing work, HiRes seed work,
    /// each noise estimate, and boundary/overlap/integration sample visits.
    pub max_work: usize,
}
impl Default for PeakPickerChromatogram {
    fn default() -> Self {
        Self {
            method: Default::default(),
            smoothing: Default::default(),
            inactive_smoothing: ChromatogramSmoothing::SavitzkyGolay {
                frame_length: 15,
                polynomial_order: 3,
            },
            peak_width: None,
            signal_to_noise: 1.0,
            // Source constructor: `updateMembers_` hands `sn_win_len` (1000)
            // and `write_sn_log_messages` (false) to `snt_`, whose own
            // defaults are 200 and true.
            noise_estimator: SignalToNoiseEstimatorMedian {
                window_length: 1000.0,
                write_log_messages: false,
                ..Default::default()
            },
            seed_signal_to_noise: 1.0,
            seed_noise_estimator: Default::default(),
            report_sn: false,
            remove_overlapping_peaks: false,
            compatibility: PickingCompatibility::default(),
            max_points: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

impl PeakPickerChromatogram {
    /// Index of the `FWHM` array in a picked chromatogram (source
    /// `FLOATINDICES::IDX_FWHM`).
    pub const IDX_FWHM: usize = 0;
    /// Index of the `IntegratedIntensity` array (source `IDX_ABUNDANCE`).
    pub const IDX_ABUNDANCE: usize = 1;
    /// Index of the `leftWidth` array (source `IDX_LEFTBORDER`).
    pub const IDX_LEFTBORDER: usize = 2;
    /// Index of the `rightWidth` array (source `IDX_RIGHTBORDER`).
    pub const IDX_RIGHTBORDER: usize = 3;
    /// Index of the `SN` array (source `IDX_SN`).
    pub const IDX_SN: usize = 4;
    /// Number of float data arrays in a picked chromatogram (source
    /// `SIZE_OF_FLOATINDICES`).
    pub const SIZE_OF_FLOATINDICES: usize = 5;

    /// The source defaults (`getDefaults`): the eleven parameters the source
    /// constructor declares, in declaration order, with their value types,
    /// descriptions and restrictions: `sgolay_frame_length` (15),
    /// `sgolay_polynomial_order` (3), `gauss_width` (50.0), `use_gauss`
    /// (`true`), `peak_width` (-1.0), `signal_to_noise` (1.0, at least 0),
    /// `sn_win_len` (1000.0), `sn_bin_count` (30), `write_sn_log_messages`
    /// (`false`), `report_sn` (`false`), `remove_overlapping_peaks` (`false`)
    /// and `method` (`corrected`; valid `legacy`, `corrected`, `crawdad`).
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure, which the constant tree cannot
    /// trigger.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this picker in the layout of
    /// [`PeakPickerChromatogram::defaults`], carrying the current values
    /// (source `getParameters`).
    ///
    /// `peak_width` is written as `-1` when [`peak_width`](Self::peak_width)
    /// is `None`; the source treats every value that is not above zero alike.
    /// The smoothing parameters come from [`smoothing`](Self::smoothing) and
    /// [`inactive_smoothing`](Self::inactive_smoothing), whichever holds that
    /// filter, and fall back to the source defaults when neither does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a count does not fit the source's
    /// 32-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let (mut gauss_width, mut frame_length, mut polynomial_order) = (50.0, 15, 3);
        for smoothing in [self.inactive_smoothing, self.smoothing] {
            match smoothing {
                ChromatogramSmoothing::Gaussian { width } => gauss_width = width,
                ChromatogramSmoothing::SavitzkyGolay {
                    frame_length: frame,
                    polynomial_order: order,
                } => {
                    frame_length = frame;
                    polynomial_order = order;
                }
            }
        }
        let int = |value: usize, key: &str| -> Result<i64> {
            i32::try_from(value)
                .map(i64::from)
                .map_err(|_| invalid(&format!("'{key}' does not fit the source 32-bit parameter")))
        };
        let entries = vec![
            int_entry(
                "sgolay_frame_length",
                int(frame_length, "sgolay_frame_length")?,
                "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.",
                false,
                None,
                None,
            ),
            int_entry(
                "sgolay_polynomial_order",
                int(polynomial_order, "sgolay_polynomial_order")?,
                "Order of the polynomial that is fitted.",
                false,
                None,
                None,
            ),
            float_entry(
                "gauss_width",
                gauss_width,
                "Gaussian width in seconds, estimated peak size.",
                false,
                None,
                None,
            ),
            string_entry(
                "use_gauss",
                if matches!(self.smoothing, ChromatogramSmoothing::Gaussian { .. }) {
                    "true"
                } else {
                    "false"
                },
                "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)",
                false,
                &["false", "true"],
            ),
            float_entry(
                "peak_width",
                self.peak_width.unwrap_or(-1.0),
                "Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.",
                false,
                None,
                None,
            ),
            float_entry(
                "signal_to_noise",
                self.signal_to_noise,
                "Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.",
                false,
                Some(0.0),
                None,
            ),
            float_entry(
                "sn_win_len",
                self.noise_estimator.window_length,
                "Signal to noise window length.",
                false,
                None,
                None,
            ),
            int_entry(
                "sn_bin_count",
                int(self.noise_estimator.bin_count, "sn_bin_count")?,
                "Signal to noise bin count.",
                false,
                None,
                None,
            ),
            flag_entry(
                "write_sn_log_messages",
                self.noise_estimator.write_log_messages,
                "Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin",
                false,
            ),
            flag_entry(
                "report_sn",
                self.report_sn,
                "Also report the apex signal-to-noise ratio for each picked peak (as an additional 'SN' FloatDataArray on the picked chromatogram). Forces initialization of the S/N estimator even if 'signal_to_noise' is 0.",
                false,
            ),
            string_entry(
                "remove_overlapping_peaks",
                if self.remove_overlapping_peaks {
                    "true"
                } else {
                    "false"
                },
                "Try to remove overlapping peaks during peak picking",
                false,
                &["false", "true"],
            ),
            string_entry(
                "method",
                self.method.name(),
                "Which method to choose for chromatographic peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed chromatogram or Crawdad on smoothed chromatogram).",
                false,
                &["legacy", "corrected", "crawdad"],
            ),
        ];
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries,
            nodes: Vec::new(),
        })
    }

    /// Build a picker from parameters, as source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing parameters take their defaults. Unknown parameters are ignored,
    /// as the source only warns; see
    /// [`PeakPickerChromatogram::from_param_with_warnings`]. The native fields
    /// ([`seed_signal_to_noise`](Self::seed_signal_to_noise),
    /// [`seed_noise_estimator`](Self::seed_noise_estimator),
    /// [`compatibility`](Self::compatibility) and the limits) keep their
    /// defaults: the source configures its seed picker once, in the
    /// constructor, and `updateMembers_` never reaches it.
    ///
    /// As in `updateMembers_`, the Savitzky-Golay and Gaussian filters are both
    /// configured whichever `use_gauss` selects, so an invalid setting of the
    /// unused filter is still refused, and `sn_win_len`, `sn_bin_count` and
    /// `write_sn_log_messages` are handed to the noise estimator's own
    /// parameters, whose restrictions (`win_len` at least one, `bin_count` at
    /// least three) then apply. A `peak_width` that is not above zero,
    /// including NaN, disables the forced width, as the source's
    /// `peak_width_ > 0.0` test does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source `Exception::InvalidParameter`);
    /// when `method` is `crawdad`, with the source's text, since this port, like
    /// the default source build, has no Crawdad backend; when the polynomial
    /// order is not below the (odd) frame length, as the Savitzky-Golay filter's
    /// `Exception::InvalidValue`; and, as native refusals, when a count is
    /// negative (the source casts it to `UInt`) or a filter setting exceeds the
    /// native limits of [`SavitzkyGolayFilter::new`] or [`GaussFilter::new`].
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`PeakPickerChromatogram::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerChromatogram::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(PEAK_PICKER_CHROMATOGRAM_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    /// Typed members from a parameter tree that already holds every default
    /// (source `updateMembers_`).
    fn from_complete_param(param: &Param) -> Result<Self> {
        let count = |key: &str| -> Result<usize> {
            usize::try_from(integer(param, key)?)
                .map_err(|_| invalid(&format!("parameter '{key}' must not be negative")))
        };
        let frame_length = count("sgolay_frame_length")?;
        let polynomial_order = count("sgolay_polynomial_order")?;
        let gauss_width = float(param, "gauss_width")?;
        let method = match param.value("method")?.as_str()? {
            "legacy" => ChromatogramPickingMethod::Legacy,
            "corrected" => ChromatogramPickingMethod::Corrected,
            "crawdad" => return Err(invalid(CRAWDAD_UNAVAILABLE)),
            _ => {
                return Err(invalid(
                    "Method needs to be one of: crawdad, corrected, legacy",
                ));
            }
        };
        // Source `sgolay_.setParameters` and `gauss_.setParameters` run for
        // both filters (`PeakPickerChromatogram.cpp:399-406`).
        SavitzkyGolayFilter::new(frame_length, polynomial_order)?;
        GaussFilter::new(GaussianWidth::Absolute(gauss_width))?;
        let gaussian = ChromatogramSmoothing::Gaussian { width: gauss_width };
        let savitzky_golay = ChromatogramSmoothing::SavitzkyGolay {
            frame_length,
            polynomial_order,
        };
        let (smoothing, inactive_smoothing) = if flag(param, "use_gauss")? {
            (gaussian, savitzky_golay)
        } else {
            (savitzky_golay, gaussian)
        };
        // Source `snt_.getParameters()`, three `setValue`s and
        // `snt_.setParameters` (`PeakPickerChromatogram.cpp:408-412`).
        let mut noise = Param::new();
        noise.set_value(
            "win_len",
            ParamValue::Float(float(param, "sn_win_len")?),
            "",
            &[],
        )?;
        noise.set_value(
            "bin_count",
            ParamValue::Integer(integer(param, "sn_bin_count")?),
            "",
            &[],
        )?;
        noise.set_value(
            "write_log_messages",
            param.value("write_sn_log_messages")?.clone(),
            "",
            &[],
        )?;
        let defaults = Self::default();
        let noise_estimator = SignalToNoiseEstimatorMedian {
            max_points: defaults.noise_estimator.max_points,
            max_bins: defaults.noise_estimator.max_bins,
            max_work: defaults.noise_estimator.max_work,
            ..SignalToNoiseEstimatorMedian::from_param(&noise)?
        };
        let peak_width = float(param, "peak_width")?;
        Ok(Self {
            method,
            smoothing,
            inactive_smoothing,
            peak_width: (peak_width > 0.0).then_some(peak_width),
            signal_to_noise: float(param, "signal_to_noise")?,
            noise_estimator,
            report_sn: flag(param, "report_sn")?,
            remove_overlapping_peaks: flag(param, "remove_overlapping_peaks")?,
            ..defaults
        })
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn spend(remaining: &mut usize, amount: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(amount)
        .ok_or_else(|| invalid("chromatogram picking work limit exceeded"))?;
    Ok(())
}

impl PeakPickerChromatogram {
    /// Pick without changing input. Five arrays are returned in source order:
    /// FWHM, IntegratedIntensity, leftWidth, rightWidth, SN. IntegratedIntensity
    /// is the inclusive sum of raw f32 samples promoted to f64, not trapezoidal area.
    ///
    /// The boundary signal is the input under
    /// [`ChromatogramPickingMethod::Legacy`] and the smoothed trace under
    /// [`ChromatogramPickingMethod::Corrected`], matching source
    /// `pickChromatogram_(chromatogram | smoothed_chrom, ...)`. Its noise is
    /// estimated with [`PickingCompatibility::source`] whatever
    /// [`compatibility`](Self::compatibility) says, because `snt_` is a plain
    /// source estimator over a signal that is either picker-generated or
    /// already validated here. A `win_len` of NaN or `+inf`, which the source's
    /// `setMinFloat("win_len", 1.0)` restriction lets through, therefore picks
    /// rather than failing.
    ///
    /// # Errors
    ///
    /// * [`Error::UnsortedData`] when retention times decrease, as source
    ///   `pickChromatogram` throws `Exception::IllegalArgument` for a
    ///   chromatogram that is not sorted by position. No compatibility flag
    ///   lifts this.
    /// * [`Error::InvalidValue`] for invalid options or resource limits,
    ///   non-finite coordinates or intensities, and — unless
    ///   [`compatibility`](Self::compatibility) allows them — negative
    ///   intensities or duplicate retention times. Also for a `win_len` below
    ///   the source's minimum of one, a `bin_count` below three, an integration
    ///   bound outside the input, an intensity that leaves `f32` range, and an
    ///   exhausted work or point budget.
    /// * The errors of the configured smoother and of the seed
    ///   [`PeakPickerHiRes`].
    pub fn pick_chromatogram(&self, input: &MSChromatogram) -> Result<ChromatogramPickingResult> {
        if self.max_points == 0
            || self.max_work == 0
            || input.len() > self.max_points
            || [self.signal_to_noise, self.seed_signal_to_noise]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || self.peak_width.is_some_and(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(invalid(
                "invalid chromatogram picker options or point limit",
            ));
        }
        input.validate()?;
        if !self.compatibility.allow_negative_intensities
            && input.peaks.iter().any(|p| p.intensity < 0.0)
        {
            return Err(invalid(
                "chromatogram picking requires nonnegative intensities; negative intensities need PickingCompatibility::allow_negative_intensities",
            ));
        }
        // The source throws `Exception::IllegalArgument` here
        // (`PeakPickerChromatogram.cpp:68-72`), so this refusal is not one
        // `allow_unsorted_positions` may lift.
        if input.peaks.windows(2).any(|p| p[0].rt > p[1].rt) {
            return Err(Error::UnsortedData);
        }
        // `isSorted` accepts equal positions, so the source reaches both the
        // seed picker and `snt_.init` with duplicate retention times.
        if !self.compatibility.allow_duplicate_positions
            && input.peaks.windows(2).any(|p| p[0].rt == p[1].rt)
        {
            return Err(invalid(
                "chromatogram picking requires distinct RT samples; duplicates need PickingCompatibility::allow_duplicate_positions",
            ));
        }
        let mut copies = super::AcquisitionCopies::default();
        copies.chromatogram(input)?;
        let mut smoothed = input.clone();
        self.smooth(&mut smoothed)?;
        let mut seed_noise = self.seed_noise_estimator.clone();
        seed_noise.max_points = seed_noise.max_points.min(self.max_points);
        seed_noise.max_work = seed_noise.max_work.min(self.max_work);
        let mut picked = PeakPickerHiRes {
            signal_to_noise: self.seed_signal_to_noise,
            noise_estimator: seed_noise,
            spacing_difference: 0.0,
            spacing_difference_gap: 0.0,
            report_fwhm: Some(FwhmUnit::Absolute),
            compatibility: self.compatibility,
            max_points: self.max_points,
            max_work: self.max_work,
            ..Default::default()
        }
        .pick_chromatogram_with_acquisition(&smoothed, false, &mut copies)?;
        let boundary_signal = match self.method {
            ChromatogramPickingMethod::Legacy => input,
            ChromatogramPickingMethod::Corrected => &smoothed,
        };
        let (noise, noise_log) = if self.signal_to_noise > 0.0 || self.report_sn {
            let mut estimator = self.noise_estimator.clone();
            estimator.max_points = estimator.max_points.min(self.max_points);
            estimator.max_work = estimator.max_work.min(self.max_work);
            // Source `snt_.init(chromatogram)`
            // (`ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:171`), where
            // `chromatogram` is this same boundary signal. `snt_` is a plain
            // source estimator with no native variant, and the signal it reads
            // is not caller data: under `corrected` this picker produced it by
            // smoothing, and under `legacy` the caller's chromatogram already
            // passed the checks above. So the estimate reproduces the source
            // unconditionally rather than following `self.compatibility`.
            //
            // Under `legacy` that is what makes a baseline-subtracted
            // chromatogram work at all: its negative samples go straight into
            // the estimate, which the strict profile refused. (Neither smoother
            // produces them from a nonnegative input: `SavitzkyGolayFilter.h:115`
            // writes `std::max(0.0, help)` and the Gaussian kernel is
            // nonnegative, so under `corrected` a negative sample survives only
            // where the local weighted average is itself negative.) It also
            // makes a NaN or infinite `win_len` pick, as the Release build
            // does, and bins an out-of-`int` histogram quotient as the Release
            // build bins it; `noise_estimator` above records why only a
            // hand-set histogram range can produce one.
            // Passing the peaks widens each `f32` with `cvtss2sd` semantics
            // rather than leaving a `f64::from` NaN payload to the compiler.
            let estimates = estimator.estimate_peaks(
                &boundary_signal.peaks,
                &PickingCompatibility::source(),
                None,
            )?;
            (Some(estimates.signal_to_noise), estimates.log)
        } else {
            (None, Vec::new())
        };
        let mut remaining = self.max_work;
        let (mut regions, apex_sn) = self.extend_boundaries(
            &boundary_signal.peaks,
            &picked.chromatogram.peaks,
            noise.as_deref(),
            &mut remaining,
        )?;
        if self.remove_overlapping_peaks {
            adjust_overlaps(&boundary_signal.peaks, &mut regions, &mut remaining)?;
        }
        let (mut areas, mut lefts, mut rights) = (Vec::new(), Vec::new(), Vec::new());
        for region in &regions {
            if region.left_index > region.right_index || region.right_index >= input.len() {
                return Err(invalid(
                    "chromatogram picker produced invalid integration bounds",
                ));
            }
            let mut area = 0.0;
            for sample in &input.peaks[region.left_index..=region.right_index] {
                spend(&mut remaining, 1)?;
                area += f64::from(sample.intensity);
            }
            areas.push(checked_intensity(area)?);
            lefts.push(checked_intensity(input.peaks[region.left_index].rt)?);
            rights.push(checked_intensity(input.peaks[region.right_index].rt)?);
        }
        // HiRes has already reported every discarded profile array. These
        // output arrays describe peaks and cannot reuse profile sample values.
        picked.chromatogram.float_data_arrays.extend([
            DataArray::new("IntegratedIntensity", areas),
            DataArray::new("leftWidth", lefts),
            DataArray::new("rightWidth", rights),
            DataArray::new("SN", apex_sn),
        ]);
        picked.chromatogram.validate()?;
        Ok(ChromatogramPickingResult {
            picked,
            smoothed,
            regions,
            noise_log,
        })
    }

    /// Replace the input only after the complete result succeeds. Use the
    /// returning method above to retain smoothing, bounds and omitted-array reports.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerChromatogram::pick_chromatogram`]; the input is then
    /// unchanged.
    pub fn filter_chromatogram(&self, input: &mut MSChromatogram) -> Result<()> {
        *input = self.pick_chromatogram(input)?.picked.chromatogram;
        Ok(())
    }

    fn smooth(&self, input: &mut MSChromatogram) -> Result<()> {
        let mut remaining = self.max_work;
        spend(&mut remaining, input.len())?;
        match self.smoothing {
            ChromatogramSmoothing::Gaussian { width } => {
                let smoother = GaussFilter::new(GaussianWidth::Absolute(width))?;
                let coefficients = (4.0 * (width / 8.0) / 0.01).ceil();
                if !coefficients.is_finite()
                    || coefficients >= smoother.algorithm.max_coefficients as f64
                {
                    return Err(invalid(
                        "Gaussian chromatogram kernel exceeds coefficient limit",
                    ));
                }
                spend(&mut remaining, coefficients as usize + 1)?;
                // The fixed 0.01 s table extends by less than 0.02 s beyond
                // width/2. Count a conservative support bound before convolution.
                let support = width / 2.0 + 0.02;
                for peak in &input.peaks {
                    let low = peak.rt - support;
                    let high = peak.rt + support;
                    if !low.is_finite() || !high.is_finite() {
                        return Err(invalid("Gaussian chromatogram support overflows"));
                    }
                    let begin = input.peaks.partition_point(|p| p.rt < low);
                    let end = input.peaks.partition_point(|p| p.rt <= high);
                    spend(&mut remaining, end - begin)?;
                }
                smoother.filter_chromatogram(input)
            }
            ChromatogramSmoothing::SavitzkyGolay {
                frame_length,
                polynomial_order,
            } => {
                let frame = frame_length
                    .checked_add(usize::from(frame_length % 2 == 0))
                    .ok_or_else(|| invalid("Savitzky-Golay frame overflows"))?;
                let setup = frame
                    .checked_mul(frame)
                    .and_then(|n| n.checked_mul(polynomial_order.checked_add(1)?))
                    .ok_or_else(|| invalid("Savitzky-Golay setup work overflows"))?;
                spend(&mut remaining, setup)?;
                if frame <= input.len() {
                    spend(
                        &mut remaining,
                        input
                            .len()
                            .checked_mul(frame)
                            .ok_or_else(|| invalid("Savitzky-Golay convolution work overflows"))?,
                    )?;
                }
                SavitzkyGolayFilter::new(frame_length, polynomial_order)?.filter_chromatogram(input)
            }
        }
    }

    fn extend_boundaries(
        &self,
        samples: &[ChromatogramPeak],
        seeds: &[ChromatogramPeak],
        sn: Option<&[f64]>,
        remaining: &mut usize,
    ) -> Result<(Vec<ChromatogramPeakRegion>, Vec<f32>)> {
        let mut regions = Vec::new();
        let mut apex_sn = Vec::new();
        let mut current = 0;
        for seed in seeds {
            current = closest_peak(samples, seed.rt, current, remaining)?;
            if current == 0 || current >= samples.len().saturating_sub(1) {
                return Err(invalid(
                    "picked apex has no two neighboring chromatogram samples",
                ));
            }
            apex_sn.push(checked_intensity(
                sn.map_or(-1.0, |values| values[current]),
            )?);
            // Source always includes the immediate neighbors. The monotonic,
            // forced-width and S/N checks begin at the second sample away.
            let (mut left, mut right) = (current - 1, current + 1);
            while left > 0 {
                spend(remaining, 1)?;
                let candidate = left - 1;
                let forced = self
                    .peak_width
                    .is_some_and(|width| (samples[candidate].rt - seed.rt).abs() < width);
                if !(samples[candidate].intensity < samples[left].intensity || forced)
                    || (self.signal_to_noise > 0.0
                        && sn.expect("positive threshold initialized noise")[candidate]
                            < self.signal_to_noise)
                {
                    break;
                }
                left = candidate;
            }
            while right + 1 < samples.len() {
                spend(remaining, 1)?;
                let candidate = right + 1;
                let forced = self
                    .peak_width
                    .is_some_and(|width| (samples[candidate].rt - seed.rt).abs() < width);
                if !(samples[candidate].intensity < samples[right].intensity || forced)
                    || (self.signal_to_noise > 0.0
                        && sn.expect("positive threshold initialized noise")[candidate]
                            < self.signal_to_noise)
                {
                    break;
                }
                right = candidate;
            }
            regions.push(ChromatogramPeakRegion {
                apex_index: current,
                left_index: left,
                right_index: right,
            });
        }
        Ok((regions, apex_sn))
    }
}

/// Source lookup intentionally selects the right neighbor on ties and returns
/// len for a target at/beyond the final sample. Callers check this end sentinel.
fn closest_peak(
    samples: &[ChromatogramPeak],
    target: f64,
    mut current: usize,
    remaining: &mut usize,
) -> Result<usize> {
    while current < samples.len() {
        spend(remaining, 1)?;
        if target < samples[current].rt {
            if current > 0
                && (target - samples[current - 1].rt).abs() < (target - samples[current].rt).abs()
            {
                current -= 1;
            }
            return Ok(current);
        }
        current += 1;
    }
    Ok(current)
}

fn adjust_overlaps(
    samples: &[ChromatogramPeak],
    regions: &mut [ChromatogramPeakRegion],
    remaining: &mut usize,
) -> Result<()> {
    for index in 0..regions.len().saturating_sub(1) {
        spend(remaining, 1)?;
        if regions[index].right_index <= regions[index + 1].left_index {
            continue;
        }
        let mut right = regions[index].apex_index;
        while right + 1 < samples.len() {
            spend(remaining, 1)?;
            if samples[right + 1].intensity >= samples[right].intensity {
                break;
            }
            right += 1;
        }
        let mut left = regions[index + 1].apex_index;
        while left > 0 {
            spend(remaining, 1)?;
            if samples[left - 1].intensity >= samples[left].intensity {
                break;
            }
            left -= 1;
        }
        if left < right {
            // These assignments are sequential in C++. The second uses the
            // updated left value and can leave overlap; do not collapse both
            // borders to a common midpoint or silently claim disjoint regions.
            left += (right - left) / 2;
            right = left + (right - left) / 2;
        }
        regions[index].right_index = right;
        regions[index + 1].left_index = left;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_closest_lookup_ties_go_right_and_final_sample_returns_end() {
        let samples = [
            ChromatogramPeak::new(1.0, 1.0),
            ChromatogramPeak::new(3.0, 1.0),
        ];
        let mut work = 100;
        assert_eq!(closest_peak(&samples, 2.0, 0, &mut work).unwrap(), 1);
        assert_eq!(closest_peak(&samples, 1.0, 0, &mut work).unwrap(), 0);
        assert_eq!(closest_peak(&samples, 0.0, 0, &mut work).unwrap(), 0);
        assert_eq!(closest_peak(&samples, 3.0, 0, &mut work).unwrap(), 2);
        assert_eq!(closest_peak(&samples, 4.0, 0, &mut work).unwrap(), 2);
        assert!(closest_peak(&samples, 2.0, 0, &mut 0).is_err());
    }

    #[test]
    fn overlap_midpoints_retain_source_sequential_assignment() {
        // Deliberately supplied seed indices model the legacy mismatch between
        // smoothed apices and a monotonically descending raw signal.
        let samples: Vec<_> = [9., 8., 7., 6., 5., 4., 3., 2., 1.]
            .iter()
            .enumerate()
            .map(|(i, &intensity)| ChromatogramPeak::new(i as f64, intensity))
            .collect();
        let mut regions = [
            ChromatogramPeakRegion {
                apex_index: 2,
                left_index: 0,
                right_index: 8,
            },
            ChromatogramPeakRegion {
                apex_index: 4,
                left_index: 3,
                right_index: 8,
            },
        ];
        adjust_overlaps(&samples, &mut regions, &mut 100).unwrap();
        assert_eq!(regions[1].left_index, 6); // (4+8)/2
        assert_eq!(regions[0].right_index, 7); // (new left6+8)/2
    }
}
