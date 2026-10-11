// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Total Ion Current (TIC) as a QC metric (source `QC/TIC.h`). See
//! `docs/QC_SUPPORT.md`.

use super::qc_base::{QCBase, Requires, Status};
use crate::Result;
use crate::format::mztab::{MzTabMetaData, MzTabParameter};
use crate::kernel::MSExperiment;
use crate::param::value::format_float32;
use crate::processing::LinearResamplerAlign;

/// Computed Total Ion Current trace plus stability counters.
///
/// Source `TIC::Result`: filled by [`Tic::compute`] and otherwise
/// default-constructed, and empty when the experiment had no spectra of the
/// requested MS level. Equality compares every field, source `operator==`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TicResult {
    /// Per-bin TIC intensities, truncated to integers.
    ///
    /// The source converts each `float` TIC to its `UInt` element, which is
    /// undefined for a value beyond `u32::MAX`; here such a value saturates
    /// at `u32::MAX` (and a negative one at 0).
    pub intensities: Vec<u32>,
    /// Per-bin intensities rescaled to `[0, 100]` of the maximum bin, and 0
    /// throughout when the maximum is 0.
    pub relative_intensities: Vec<f32>,
    /// Per-bin retention times in seconds, narrowed to `f32` as the source
    /// stores them.
    pub retention_times: Vec<f32>,
    /// Sum of [`TicResult::intensities`], the area under the TIC.
    ///
    /// The source sums into a `UInt`, so an area beyond `u32::MAX`
    /// (4,294,967,295) wraps modulo 2^32 there. This port sums exactly into
    /// a `u64`; the two agree whenever the source does not wrap.
    pub area: u64,
    /// Number of consecutive bins where the intensity drops by at least a
    /// factor of ten: `intensities[i] < intensities[i - 1] / 10`, with integer
    /// division.
    pub fall: u32,
    /// Number of consecutive bins where the intensity rises by at least a
    /// factor of ten: `intensities[i] > intensities[i - 1] * 10`.
    ///
    /// The source multiplies in `UInt`, which wraps once the earlier bin
    /// exceeds 429,496,729 and then miscounts; this port multiplies exactly.
    pub jump: u32,
}

/// Total Ion Current (TIC) as a QC metric.
///
/// Source `TIC`: builds a TIC trace at a chosen MS level and returns a
/// [`TicResult`] with the per-bin intensities, the relative intensities
/// (0 to 100), the retention times, the area under the TIC and two stability
/// counters for abrupt ten-fold changes between consecutive bins. Each call to
/// [`Tic::compute`] is independent; nothing is retained by the instance.
///
/// The source header also declares `getResults()`, returning chromatograms,
/// but never defines it, so a call does not link; it is not ported.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tic;

impl Tic {
    /// Source `getName`: `TIC`.
    pub const NAME: &'static str = "TIC";

    /// A new metric, source `TIC()`.
    pub fn new() -> Self {
        Self
    }

    /// Compute the Total Ion Current trace and stability metrics, source
    /// `compute(exp, bin_size = 0, ms_level = 1)`.
    ///
    /// The trace is `MSExperiment::calculateTIC(bin_size, ms_level)`: one
    /// point per spectrum of MS level `ms_level` (every level when it is 0), in
    /// spectrum order, at the spectrum's retention time with its summed
    /// intensity. When `bin_size` is positive, the trace is then resampled
    /// onto a grid of `bin_size` seconds with [`LinearResamplerAlign`]; 0 (or
    /// any non-positive or NaN value) disables binning.
    ///
    /// The result is empty when `exp` has no spectrum of MS level `ms_level`.
    ///
    /// # Errors
    ///
    /// Only with binning: as [`LinearResamplerAlign::new`] for a non-finite
    /// `bin_size` and [`LinearResamplerAlign::raster_chromatogram`] for a
    /// trace it cannot resample. The source has no error path.
    pub fn compute(&self, exp: &MSExperiment, bin_size: f32, ms_level: u32) -> Result<TicResult> {
        let mut trace = exp.calculate_tic(ms_level);
        if bin_size > 0.0 {
            LinearResamplerAlign::new(f64::from(bin_size))?.raster_chromatogram(&mut trace)?;
        }
        let mut result = TicResult::default();
        if trace.peaks.is_empty() {
            return Ok(result);
        }
        for peak in &trace.peaks {
            // `as` saturates and maps NaN to 0: the documented narrowing.
            result.intensities.push(peak.intensity as u32);
            result.retention_times.push(peak.rt as f32);
        }
        let max_int = result.intensities.iter().copied().max().unwrap_or(0);
        result.relative_intensities = result
            .intensities
            .iter()
            .map(|&i| {
                if max_int != 0 {
                    (f64::from(i) / f64::from(max_int) * 100.0) as f32
                } else {
                    0.0
                }
            })
            .collect();
        result.area = u64::from(result.intensities[0]);
        for pair in result.intensities.windows(2) {
            let (previous, current) = (u64::from(pair[0]), u64::from(pair[1]));
            result.area += current;
            if current > previous * 10 {
                result.jump += 1;
            }
            if current < previous / 10 {
                result.fall += 1;
            }
        }
        Ok(result)
    }

    /// Append the TIC traces of `tics` to the metadata section `meta` as
    /// custom parameters, source `addMetaDataMetricsToMzTab`.
    ///
    /// For every non-empty entry of `tics`, one parameter is added with the
    /// controlled-vocabulary label `total ion current` (sic: the source puts
    /// the term's name where the vocabulary label belongs), accession
    /// `MS:1000285`, the name `TIC_<n>` with `n` the 1-based index into
    /// `tics`, and the value `[rt0, int0, rt1, int1, ...]`: each retention time
    /// written as the source's `StringUtils::toStr(float, false)` writes it
    /// (three fractional digits, trailing zeros trimmed) and each intensity as
    /// an integer. Empty entries (no spectra of that level) are skipped.
    ///
    /// Each parameter is stored under the key `meta.custom.len()`, the
    /// source's `meta.custom[meta.custom.size()]`: 0-based, and when the map's
    /// keys are not already `0..len` the new entry replaces the one at that
    /// key, as in the source.
    pub fn add_meta_data_metrics_to_mz_tab(&self, meta: &mut MzTabMetaData, tics: &[TicResult]) {
        for (index, tic) in tics.iter().enumerate() {
            if tic.intensities.is_empty() {
                continue; // no spectra of this level
            }
            let mut value = String::from("[");
            for (position, (rt, intensity)) in
                tic.retention_times.iter().zip(&tic.intensities).enumerate()
            {
                if position > 0 {
                    value.push_str(", ");
                }
                value.push_str(&format_float32(*rt, false));
                value.push_str(", ");
                value.push_str(&intensity.to_string());
            }
            value.push(']');
            let parameter = MzTabParameter::from_parts(
                "total ion current",
                "MS:1000285",
                format!("TIC_{}", index + 1),
                value,
            );
            let key = meta.custom.len();
            meta.custom.insert(key, parameter);
        }
    }
}

impl QCBase for Tic {
    fn name(&self) -> &str {
        Self::NAME
    }
    /// Raw mzML data, source `requirements`.
    fn requirements(&self) -> Status {
        Status::from_flag(Requires::RawMzml)
    }
}
