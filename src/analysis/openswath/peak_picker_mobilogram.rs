// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Peak picking in a single ion mobilogram.
//!
//! Ported from `ANALYSIS/OPENSWATH/PeakPickerMobilogram.h` and `.cpp` at
//! OpenMS4-core `bc9cc12`. The source class is `PeakPickerChromatogram` on
//! ion mobility instead of retention time, without the boundary S/N gate:
//! smoothing, `PeakPickerHiRes` seeding, boundary extension, overlap removal
//! and raw-intensity integration are the same template code on a `Mobilogram`
//! (`GaussFilter::filter(Mobilogram&)`, `SavitzkyGolayFilter::filter`,
//! `PeakPickerHiRes::pick(const Mobilogram&, ...)`). The port therefore runs
//! [`PeakPickerChromatogram`] on the mobilogram's (mobility, intensity) samples
//! with the S/N gate off, and adds the mobilogram-only parts: the two index
//! arrays and the top-peak filters. See
//! `docs/MRM_TRANSITION_GROUP_PICKER_SUPPORT.md` (section "Mobilogram peak
//! picking").
//!
//! [`PeakPickerChromatogram`]: crate::processing::chromatogram::PeakPickerChromatogram

use crate::kernel::{ChromatogramPeak, DataArray, MSChromatogram, MobilityPeak1D, Mobilogram};
use crate::param::{DefaultParamHandler, Param, ParamNode};
use crate::processing::chromatogram::{
    ChromatogramPeakRegion, ChromatogramPickingMethod, ChromatogramSmoothing,
    PeakPickerChromatogram,
};
use crate::processing::peak_picking::{
    PickingCompatibility, flag, flag_entry, float, float_entry, int_entry, integer, string_entry,
};
use crate::processing::smoothing::{GaussFilter, GaussianWidth, SavitzkyGolayFilter};
use crate::{Error, Result};

/// The source class name, `DefaultParamHandler("PeakPickerMobilogram")`,
/// under which unknown-parameter warnings are reported.
pub const PEAK_PICKER_MOBILOGRAM_NAME: &str = "PeakPickerMobilogram";

/// Left, apex and right position of the highest picked peak (source
/// `PeakPickerMobilogram::PeakPositions`).
///
/// `left` and `right` are indices into the mobilogram. `apex` is what the
/// source's `findHighestPeak_` returns, which is the index of the picked peak
/// in the picked mobilogram, not an index into the mobilogram, although the
/// source documents all three as mobilogram indices; this port keeps the
/// source value.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PeakPositions {
    /// Left position of the peak.
    pub left: usize,
    /// Apex position of the peak (see the type documentation).
    pub apex: usize,
    /// Right position of the peak.
    pub right: usize,
}

/// What one [`PeakPickerMobilogram::pick_mobilogram`] call produces: the
/// source's `picked_mobilogram` and `smoothed_mobilogram` out-parameters and
/// its three public per-peak vectors.
#[derive(Clone, Debug, PartialEq)]
pub struct MobilogramPickingResult {
    /// The picked peaks with the six float arrays at the `IDX_*` indices:
    /// `FWHM`, `IntegratedIntensity`, `leftWidth`, `rightWidth`,
    /// `leftWidthIndex`, `rightWidthIndex`.
    pub picked: Mobilogram,
    /// The smoothed mobilogram.
    pub smoothed: Mobilogram,
    /// Source member `integrated_intensities_`: the inclusive sum of the raw
    /// intensities of each peak.
    pub integrated_intensities: Vec<f64>,
    /// Source member `left_width_`: the first sample of each peak.
    pub left_widths: Vec<usize>,
    /// Source member `right_width_`: the last sample of each peak.
    pub right_widths: Vec<usize>,
}

/// Finds peaks in a single mobilogram and annotates their left and right
/// borders (source `PeakPickerMobilogram`).
///
/// It uses the `PeakPickerHiRes` internally to find seed candidates, which are
/// then expanded and a left and right border searched. Additionally,
/// overlapping peaks can be removed.
///
/// The source class is a `DefaultParamHandler`; its eleven parameters are
/// [`PeakPickerMobilogram::defaults`], [`PeakPickerMobilogram::to_param`] and
/// [`PeakPickerMobilogram::from_param`].
#[derive(Clone, Debug)]
pub struct PeakPickerMobilogram {
    /// Source `method`: `legacy` (boundaries on the raw mobilogram) or
    /// `corrected` (the default; boundaries on the smoothed one). Both
    /// integrate the raw intensities.
    pub method: ChromatogramPickingMethod,
    /// Source `use_gauss` (default `false`) with `sgolay_frame_length` (9) and
    /// `sgolay_polynomial_order` (3), or `gauss_width` (0.002).
    pub smoothing: ChromatogramSmoothing,
    /// The parameters of the smoother `smoothing` does not select, kept for
    /// [`to_param`](Self::to_param); see
    /// [`PeakPickerChromatogram::inactive_smoothing`].
    pub inactive_smoothing: ChromatogramSmoothing,
    /// Source `peak_width`: extend the peak at least this far on both sides;
    /// `None` (the default, `-1`) turns this off.
    pub peak_width: Option<f64>,
    /// Source `signal_to_noise` (default 1.0, at least 0). The source stores
    /// it and uses it for nothing after construction: the boundary walk has no
    /// S/N gate, and the seed picker took the constructor's value.
    pub signal_to_noise: f64,
    /// Source `sn_win_len` (an integer parameter, default 1); stored, no
    /// effect, as in the source.
    pub sn_win_len: i64,
    /// Source `sn_bin_count` (default 4); stored, no effect, as in the source.
    pub sn_bin_count: i64,
    /// Source `write_sn_log_messages` (default `false`); stored, no effect, as
    /// in the source.
    pub write_sn_log_messages: bool,
    /// Source `remove_overlapping_peaks` (default `false`).
    ///
    /// The source declares this parameter but `updateMembers_` never assigns
    /// its member `remove_overlapping_`, which `pickMobilogram` then reads
    /// uninitialised; this port assigns it from the parameter, as the
    /// parameter's description promises.
    pub remove_overlapping_peaks: bool,
    /// S/N threshold of the seed `PeakPickerHiRes`: 1, the default
    /// `signal_to_noise` the source constructor copies into `pp_`.
    pub seed_signal_to_noise: f64,
    /// Source behaviours the native default refuses; see
    /// [`PeakPickerChromatogram::compatibility`].
    pub compatibility: PickingCompatibility,
    /// Native ceiling on the input point count.
    pub max_points: usize,
    /// Native per-stage work limit; see
    /// [`PeakPickerChromatogram::max_work`].
    pub max_work: usize,
}

impl Default for PeakPickerMobilogram {
    fn default() -> Self {
        Self {
            method: ChromatogramPickingMethod::Corrected,
            smoothing: ChromatogramSmoothing::SavitzkyGolay {
                frame_length: 9,
                polynomial_order: 3,
            },
            inactive_smoothing: ChromatogramSmoothing::Gaussian { width: 0.002 },
            peak_width: None,
            signal_to_noise: 1.0,
            sn_win_len: 1,
            sn_bin_count: 4,
            write_sn_log_messages: false,
            remove_overlapping_peaks: false,
            seed_signal_to_noise: 1.0,
            compatibility: PickingCompatibility::default(),
            max_points: 1_000_000,
            max_work: 50_000_000,
        }
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

impl PeakPickerMobilogram {
    /// Index of the `FWHM` array (source `FLOATINDICES::IDX_FWHM`).
    pub const IDX_FWHM: usize = 0;
    /// Index of the `IntegratedIntensity` array (source `IDX_ABUNDANCE`).
    pub const IDX_ABUNDANCE: usize = 1;
    /// Index of the `leftWidth` array (source `IDX_LEFTBORDER`).
    pub const IDX_LEFTBORDER: usize = 2;
    /// Index of the `rightWidth` array (source `IDX_RIGHTBORDER`).
    pub const IDX_RIGHTBORDER: usize = 3;
    /// Index of the `leftWidthIndex` array (source `IDX_OF_LEFTBORDER_IDX`).
    pub const IDX_OF_LEFTBORDER_IDX: usize = 4;
    /// Index of the `rightWidthIndex` array (source `IDX_OF_RIGHTBORDER_IDX`).
    pub const IDX_OF_RIGHTBORDER_IDX: usize = 5;
    /// Number of float arrays of a picked mobilogram (source
    /// `SIZE_OF_FLOATINDICES`).
    pub const SIZE_OF_FLOATINDICES: usize = 6;

    /// The source defaults (`getDefaults`), in declaration order:
    /// `sgolay_frame_length` (9), `sgolay_polynomial_order` (3), `gauss_width`
    /// (0.002), `use_gauss` (`false`), `peak_width` (-1.0), `signal_to_noise`
    /// (1.0, at least 0), `sn_win_len` (the integer 1), `sn_bin_count` (4),
    /// `write_sn_log_messages` (`false`), `remove_overlapping_peaks`
    /// (`false`) and `method` (`corrected`; valid `legacy`, `corrected`,
    /// `crawdad`), with the source descriptions.
    ///
    /// # Errors
    ///
    /// Only a parameter-tree resource failure.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this picker in the layout of
    /// [`PeakPickerMobilogram::defaults`], carrying the current values
    /// (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a smoothing count does not fit the
    /// source's 32-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let (mut gauss_width, mut frame_length, mut polynomial_order) = (0.002, 9, 3);
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
            int_entry(
                "sn_win_len",
                self.sn_win_len,
                "Signal to noise window length.",
                false,
                None,
                None,
            ),
            int_entry(
                "sn_bin_count",
                self.sn_bin_count,
                "Signal to noise bin count.",
                false,
                None,
                None,
            ),
            flag_entry(
                "write_sn_log_messages",
                self.write_sn_log_messages,
                "Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin",
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
                "Which method to choose for mobilogram peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed mobilogram).",
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
    /// Missing parameters take their defaults; unknown ones are ignored, as
    /// the source only warns (see
    /// [`PeakPickerMobilogram::from_param_with_warnings`]). Both smoothers are
    /// configured whichever `use_gauss` selects, as `updateMembers_` does. The
    /// native fields keep their defaults.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a wrong value type or a violated
    /// restriction (the source `Exception::InvalidParameter`); for a
    /// Savitzky-Golay order not below the frame length (the filter's
    /// `Exception::InvalidValue`); for a negative smoothing count (the source
    /// casts it to `UInt`); and for `method` `crawdad`, which the source
    /// accepts and then handles in neither branch of `pickMobilogram`, reading
    /// the per-peak vectors of an earlier call (or none) — a native refusal.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`PeakPickerMobilogram::from_param`], also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`PeakPickerMobilogram::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(PEAK_PICKER_MOBILOGRAM_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

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
            other => {
                return Err(invalid(&format!(
                    "PeakPickerMobilogram method '{other}' is not supported: the source handles only legacy and corrected"
                )));
            }
        };
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
        let peak_width = float(param, "peak_width")?;
        Ok(Self {
            method,
            smoothing,
            inactive_smoothing,
            peak_width: (peak_width > 0.0).then_some(peak_width),
            signal_to_noise: float(param, "signal_to_noise")?,
            sn_win_len: integer(param, "sn_win_len")?,
            sn_bin_count: integer(param, "sn_bin_count")?,
            write_sn_log_messages: flag(param, "write_sn_log_messages")?,
            remove_overlapping_peaks: flag(param, "remove_overlapping_peaks")?,
            ..Self::default()
        })
    }

    /// Find peaks in a single mobilogram and annotate their left and right
    /// borders (source `pickMobilogram`, both overloads).
    ///
    /// The mobilogram is smoothed, seeded with the `PeakPickerHiRes`, and each
    /// seed's borders are extended on the raw (`legacy`) or smoothed
    /// (`corrected`) signal while intensities keep falling or the forced
    /// [`peak_width`](Self::peak_width) reaches; overlapping neighbours are
    /// then split when [`remove_overlapping_peaks`](Self::remove_overlapping_peaks)
    /// is set, and each peak's raw intensities are summed. Unlike
    /// `PeakPickerChromatogram` there is no S/N gate and no `SN` array; the
    /// two index arrays hold the border sample indices as `f32`.
    ///
    /// An empty mobilogram returns empty results, as the source returns
    /// before touching its outputs.
    ///
    /// # Errors
    ///
    /// [`Error::UnsortedData`] when mobilities decrease (the source's
    /// `Exception::IllegalArgument`), and the errors of
    /// [`PeakPickerChromatogram::pick_chromatogram`] on the same samples.
    pub fn pick_mobilogram(&self, input: &Mobilogram) -> Result<MobilogramPickingResult> {
        if input.peaks.is_empty() {
            return Ok(MobilogramPickingResult {
                picked: Mobilogram::default(),
                smoothed: Mobilogram::default(),
                integrated_intensities: Vec::new(),
                left_widths: Vec::new(),
                right_widths: Vec::new(),
            });
        }
        if input.len() > self.max_points {
            return Err(invalid("mobilogram exceeds the peak picker point limit"));
        }
        input.validate()?;
        let trace = MSChromatogram {
            peaks: input
                .peaks
                .iter()
                .map(|p| ChromatogramPeak::new(p.mobility, p.intensity))
                .collect(),
            ..MSChromatogram::default()
        };
        let picker = PeakPickerChromatogram {
            method: self.method,
            smoothing: self.smoothing,
            inactive_smoothing: self.inactive_smoothing,
            peak_width: self.peak_width,
            signal_to_noise: 0.0,
            seed_signal_to_noise: self.seed_signal_to_noise,
            report_sn: false,
            remove_overlapping_peaks: self.remove_overlapping_peaks,
            compatibility: self.compatibility,
            max_points: self.max_points,
            max_work: self.max_work,
            ..PeakPickerChromatogram::default()
        };
        let result = picker.pick_chromatogram(&trace)?;
        let picked_trace = &result.picked.chromatogram;
        let array = |index: usize| -> Result<Vec<f32>> {
            picked_trace
                .float_data_arrays
                .get(index)
                .map(|a| a.data.clone())
                .ok_or_else(|| invalid("chromatogram picker returned too few arrays"))
        };
        let regions: &[ChromatogramPeakRegion] = &result.regions;
        let left_widths: Vec<usize> = regions.iter().map(|r| r.left_index).collect();
        let right_widths: Vec<usize> = regions.iter().map(|r| r.right_index).collect();
        let integrated_intensities: Vec<f64> = regions
            .iter()
            .map(|r| {
                input.peaks[r.left_index..=r.right_index]
                    .iter()
                    .fold(0.0, |sum, p| sum + f64::from(p.intensity))
            })
            .collect();
        let picked = Mobilogram {
            peaks: picked_trace
                .peaks
                .iter()
                .map(|p| MobilityPeak1D::new(p.rt, p.intensity))
                .collect(),
            rt: input.rt,
            drift_time_unit: input.drift_time_unit,
            float_data_arrays: vec![
                DataArray::new("FWHM", array(PeakPickerChromatogram::IDX_FWHM)?),
                DataArray::new(
                    "IntegratedIntensity",
                    array(PeakPickerChromatogram::IDX_ABUNDANCE)?,
                ),
                DataArray::new("leftWidth", array(PeakPickerChromatogram::IDX_LEFTBORDER)?),
                DataArray::new(
                    "rightWidth",
                    array(PeakPickerChromatogram::IDX_RIGHTBORDER)?,
                ),
                DataArray::new(
                    "leftWidthIndex",
                    left_widths.iter().map(|&i| i as f32).collect(),
                ),
                DataArray::new(
                    "rightWidthIndex",
                    right_widths.iter().map(|&i| i as f32).collect(),
                ),
            ],
            integer_data_arrays: Vec::new(),
            string_data_arrays: Vec::new(),
        };
        picked.validate()?;
        let smoothed = Mobilogram {
            peaks: result
                .smoothed
                .peaks
                .iter()
                .map(|p| MobilityPeak1D::new(p.rt, p.intensity))
                .collect(),
            ..input.clone()
        };
        Ok(MobilogramPickingResult {
            picked,
            smoothed,
            integrated_intensities,
            left_widths,
            right_widths,
        })
    }

    /// Keep only the samples of the highest picked peak in each of
    /// `mobilograms` (source `filterTopPeak`, vector overload), and return its
    /// positions.
    ///
    /// The highest peak is the one with the largest `IntegratedIntensity` in
    /// `picked` (the first on ties); every mobilogram is cut to its
    /// `[leftWidthIndex, rightWidthIndex]` sample range, keeping parallel
    /// arrays aligned. A mobilogram the range does not fit is left unchanged.
    /// Without picked peaks the range is the whole first mobilogram (`apex`
    /// its middle); without the six arrays, or without mobilograms, the
    /// result is all zeros and nothing is cut.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the abundance and index arrays differ in
    /// length (the source reads out of bounds), and the errors of
    /// [`Mobilogram::select`].
    pub fn filter_top_peak(
        picked: &Mobilogram,
        mobilograms: &mut [Mobilogram],
    ) -> Result<PeakPositions> {
        if picked.float_data_arrays.len() <= Self::IDX_OF_RIGHTBORDER_IDX || mobilograms.is_empty()
        {
            return Ok(PeakPositions::default());
        }
        let positions = Self::highest_peak(picked, mobilograms[0].len())?;
        for mobilogram in mobilograms.iter_mut() {
            filter_peak_intensities(mobilogram, positions.left, positions.right)?;
        }
        Ok(positions)
    }

    /// [`PeakPickerMobilogram::filter_top_peak`] for one mobilogram (source
    /// `filterTopPeak`, single overload).
    ///
    /// # Errors
    ///
    /// As [`PeakPickerMobilogram::filter_top_peak`].
    pub fn filter_top_peak_single(
        picked: &Mobilogram,
        mobilogram: &mut Mobilogram,
    ) -> Result<PeakPositions> {
        if picked.float_data_arrays.len() <= Self::IDX_OF_RIGHTBORDER_IDX {
            return Ok(PeakPositions::default());
        }
        let positions = Self::highest_peak(picked, mobilogram.len())?;
        filter_peak_intensities(mobilogram, positions.left, positions.right)?;
        Ok(positions)
    }

    /// Source `findHighestPeak_` over the picked arrays (after
    /// `extractFloatValues_` and `extractIntValues_`).
    fn highest_peak(picked: &Mobilogram, im_size: usize) -> Result<PeakPositions> {
        let abundances = &picked.float_data_arrays[Self::IDX_ABUNDANCE].data;
        let lefts = &picked.float_data_arrays[Self::IDX_OF_LEFTBORDER_IDX].data;
        let rights = &picked.float_data_arrays[Self::IDX_OF_RIGHTBORDER_IDX].data;
        if abundances.is_empty() {
            return Ok(PeakPositions {
                left: 0,
                apex: im_size / 2,
                right: im_size.wrapping_sub(1),
            });
        }
        let mut max_index = 0;
        for (i, value) in abundances.iter().enumerate() {
            if f64::from(abundances[max_index]) < f64::from(*value) {
                max_index = i;
            }
        }
        let index = |values: &[f32]| -> Result<usize> {
            values
                .get(max_index)
                .map(|&v| float_to_size(v))
                .ok_or_else(|| {
                    invalid("picked mobilogram index arrays are shorter than its abundances")
                })
        };
        Ok(PeakPositions {
            left: index(lefts)?,
            apex: max_index,
            right: index(rights)?,
        })
    }
}

/// `static_cast<std::size_t>(float)` as the Linux x86-64 build performs it for
/// the stored border indices: truncation toward zero; a negative or NaN value,
/// which the picker never stores, maps to zero here rather than to the
/// target's out-of-range result.
fn float_to_size(value: f32) -> usize {
    if value.is_nan() || value <= 0.0 {
        0
    } else {
        value as usize
    }
}

/// Source `filterPeakIntensities_`: keep the samples `left..=right`, or leave
/// the mobilogram unchanged when it is empty or the range does not fit.
fn filter_peak_intensities(mobilogram: &mut Mobilogram, left: usize, right: usize) -> Result<()> {
    if mobilogram.is_empty() || right >= mobilogram.len() || left > right {
        return Ok(());
    }
    let keep: Vec<usize> = (left..=right).collect();
    mobilogram.select(&keep)
}
