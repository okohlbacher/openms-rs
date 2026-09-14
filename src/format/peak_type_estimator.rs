// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Spectrum peak-type estimation, the port of `FORMAT/PeakTypeEstimator.h`.
//!
//! The source header is header-only: one static template,
//! `PeakTypeEstimator::estimateType(begin, end)`, which classifies the peaks of
//! an iterator range as profile (raw) or centroided (peak) data. This module is
//! its public entry point, `PeakTypeEstimator::estimate_type`, over a slice of
//! `Peak1D`. The shoulder algorithm itself is the line-by-line
//! transcription that `MSSpectrum::getType(true)` already uses
//! (`src/kernel/spectrum_type.rs`), so the spectrum query and this function
//! can never disagree on the same peaks.
//!
//! See `docs/PEAK_TYPE_ESTIMATOR_SUPPORT.md` for the API mapping, the preserved
//! source conventions, the native differences and the evidence.
//!
//! There is a second, stricter classifier in the crate,
//! `processing::peak_picking::estimate_spectrum_type`. It validates the whole
//! spectrum record and rejects negative intensities and unsorted or duplicate
//! m/z values before classifying. The source estimator checks none of that, and
//! FileInfo and `MSSpectrum::getType` call the source estimator, so source
//! callers belong here.

use crate::kernel::{Peak1D, SpectrumType, SpectrumTypeQueryLimits};
use crate::{Error, Result};

/// Work units charged per peak, the same charge `MSSpectrum::get_type_with_limits`
/// makes before it estimates: two scalar copies and finite checks, the total
/// sum and the five complete maximum and shoulder scans.
const WORK_PER_POINT: usize = 32;

/// Bytes of scratch storage per peak: one `f64` m/z copy and one `f64`
/// intensity copy.
const BYTES_PER_POINT: usize = 2 * std::mem::size_of::<f64>();

fn resource() -> Error {
    Error::InvalidValue("spectrum type query resource limit exceeded".into())
}

/// Estimates whether the data of a spectrum are raw (profile) data or peak
/// (centroided) data, the source `class PeakTypeEstimator`.
///
/// The source class holds no state and exists only to scope its one static
/// template; its class test nevertheless constructs and deletes an instance
/// (`PeakTypeEstimator_test.cpp:33-40`, both `[EXTRA]`). This is therefore a
/// zero-sized unit struct: [`Default`] is the constructor, and it has no drop
/// glue, which is the destructor.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PeakTypeEstimator;

impl PeakTypeEstimator {
    /// Fewest peaks the estimator classifies; shorter input is
    /// [`SpectrumType::Unknown`].
    ///
    /// The source `@note`: *if there are less than 5 peaks in the iterator
    /// range `SpectrumType::UNKNOWN` is returned* (`PeakTypeEstimator.h:40`,
    /// checked at line 47).
    pub const MIN_PEAKS: usize = 5;

    /// Estimates the peak type of `peaks` from the intensity characteristics of
    /// up to five maxima, the source `estimateType(begin, end)`.
    ///
    /// The source documentation, carried across: profile and centroided data
    /// are told apart by looking at the highest five peaks of the spectrum. If
    /// at least two neighbouring sampling points on either side of a local
    /// maximum lie within 1 Th, the peak counts as a profile peak; the
    /// intensities must decline on both shoulders. Every sampling point
    /// assigned to a shoulder is excluded from the search for the next highest
    /// peak, except a local minimum at the end of a shoulder, which may serve
    /// the shoulders of the peaks to its left and right. Once five peaks or
    /// more than 50% of the total spectral intensity have been examined, the
    /// number of peaks classified as centroided (C) is compared with the
    /// number classified as profile (P): the spectrum is
    /// [`SpectrumType::Profile`] if P / (C + P) > 0.75 and
    /// [`SpectrumType::Centroid`] otherwise.
    ///
    /// Fewer than [`MIN_PEAKS`](Self::MIN_PEAKS) peaks give
    /// [`SpectrumType::Unknown`] without looking at their values, as in the
    /// source.
    ///
    /// Source conventions this preserves (see the support document for each):
    ///
    /// - A shoulder point qualifies only while its intensity is positive, at
    ///   most the previous point's, *strictly* more than 10% of it, and
    ///   *strictly* less than 1 Th from the maximum.
    /// - No sortedness, distinctness or sign requirement: unsorted m/z,
    ///   duplicate m/z and negative intensities are classified, not refused.
    /// - A spectrum with no positive intensity yields no evidence at all, and
    ///   the source's `0 / 0.0f` evidence ratio is NaN, which is not greater
    ///   than 0.75: such input is [`SpectrumType::Centroid`].
    /// - The source copies the peaks and works on the copy; `peaks` is borrowed
    ///   immutably here and the scratch copy is internal.
    ///
    /// The source works on `float` intensities promoted to `double`; this works
    /// on the same values held in `f64`. Every value the algorithm writes back
    /// is either zero or an intensity it read, so the scratch values stay
    /// exactly the input `f32` values, and every comparison and the running
    /// `double` sums are the source's, operation for operation.
    ///
    /// A sub-range of a spectrum, which the source expresses with an iterator
    /// pair, is a sub-slice here: `estimate_type(&spectrum.peaks[..4])`
    /// corresponds to `estimateType(spec.begin(), spec.begin() + 4)`.
    ///
    /// Uses [`SpectrumTypeQueryLimits::default`]; see
    /// [`estimate_type_with_limits`](Self::estimate_type_with_limits) for the
    /// ceilings and the errors.
    ///
    /// # Errors
    ///
    /// As [`estimate_type_with_limits`](Self::estimate_type_with_limits) with
    /// the default limits.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::format::peak_type_estimator::PeakTypeEstimator;
    /// use openms::kernel::{Peak1D, SpectrumType};
    ///
    /// // One profile peak sampled every 0.01 Th.
    /// let profile: Vec<Peak1D> = [10.0, 30.0, 60.0, 100.0, 60.0, 30.0, 10.0]
    ///     .iter()
    ///     .enumerate()
    ///     .map(|(i, &intensity)| Peak1D::new(500.0 + 0.01 * i as f64, intensity))
    ///     .collect();
    /// assert_eq!(PeakTypeEstimator::estimate_type(&profile)?, SpectrumType::Profile);
    /// // Fewer than five points are never classified.
    /// assert_eq!(PeakTypeEstimator::estimate_type(&profile[..4])?, SpectrumType::Unknown);
    /// # Ok::<(), openms::Error>(())
    /// ```
    pub fn estimate_type(peaks: &[Peak1D]) -> Result<SpectrumType> {
        Self::estimate_type_with_limits(peaks, SpectrumTypeQueryLimits::default())
    }

    /// As [`estimate_type`](Self::estimate_type), with explicit resource
    /// ceilings.
    ///
    /// The ceilings and their order are exactly those
    /// `MSSpectrum::get_type_with_limits` applies to the peaks once it reaches
    /// estimation:
    ///
    /// 1. fewer than [`MIN_PEAKS`](Self::MIN_PEAKS) peaks return
    ///    [`SpectrumType::Unknown`] before any ceiling is consulted;
    /// 2. more than `limits.max_points` peaks fail;
    /// 3. 32 work units per peak must fit in `limits.max_work`;
    /// 4. 16 bytes of scratch storage per peak must fit in `limits.max_bytes`;
    /// 5. every m/z and intensity must be finite;
    /// 6. only then is the scratch copy allocated, fallibly.
    ///
    /// Queried with `query_data` set, a spectrum with no stored type and no
    /// data-processing records therefore succeeds or fails together with its
    /// peak slice under the same limits, with the same class. A spectrum that
    /// carries data-processing records need not: while `get_type_with_limits`
    /// searches the records for a peak-picking step, it charges each one
    /// `1 + 12 * h` work units against the same `max_work` before the peaks
    /// are charged, where `h` is the bit length of the record's action count
    /// (one unit for a record without actions). Such a spectrum can exceed
    /// `max_work` under limits its peak slice fits: seven peaks at `max_work`
    /// 224 are classified as a slice and refused as a spectrum carrying one
    /// default record. A record with a peak-picking action instead makes the
    /// spectrum [`SpectrumType::Centroid`] without estimating. This function
    /// sees no records and charges none.
    ///
    /// The whole cost is linear in the number of peaks: at most five maxima are
    /// examined and each shoulder scan clears the points it passes.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a ceiling is exceeded, when the
    /// scratch allocation fails, and when any m/z or intensity is NaN or
    /// infinite. The source checks none of these; non-finite values would make
    /// its comparisons silently meaningless (a NaN intensity is never a
    /// maximum and poisons the total), so this refuses them instead of
    /// returning a class derived from them. Nothing is mutated in any case.
    pub fn estimate_type_with_limits(
        peaks: &[Peak1D],
        limits: SpectrumTypeQueryLimits,
    ) -> Result<SpectrumType> {
        let n = peaks.len();
        if n < Self::MIN_PEAKS {
            return Ok(SpectrumType::Unknown);
        }
        if n > limits.max_points {
            return Err(resource());
        }
        let work = n.checked_mul(WORK_PER_POINT).ok_or_else(resource)?;
        if work > limits.max_work {
            return Err(resource());
        }
        let bytes = n.checked_mul(BYTES_PER_POINT).ok_or_else(resource)?;
        if bytes > limits.max_bytes {
            return Err(resource());
        }
        if peaks
            .iter()
            .any(|peak| !peak.mz.is_finite() || !peak.intensity.is_finite())
        {
            return Err(Error::InvalidValue(
                "spectrum type estimation requires finite consumed peak values".into(),
            ));
        }
        let mut mz = Vec::new();
        let mut intensity = Vec::new();
        mz.try_reserve_exact(n).map_err(|_| resource())?;
        intensity.try_reserve_exact(n).map_err(|_| resource())?;
        for peak in peaks {
            mz.push(peak.mz);
            intensity.push(f64::from(peak.intensity));
        }
        Ok(crate::kernel::spectrum_type::estimate(&mz, &mut intensity))
    }
}
