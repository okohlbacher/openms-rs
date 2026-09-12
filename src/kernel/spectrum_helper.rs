// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Free helper functions for spectra and chromatograms.
//!
//! Source: `KERNEL/SpectrumHelper.h` and `KERNEL/SpectrumHelper.cpp` at Core SDK
//! `bc9cc12`. The support document `docs/SPECTRUM_HELPER_SUPPORT.md` lists every
//! source member, its counterpart here, the preserved source conventions and the
//! native differences. Source `removePeaks` is not duplicated; it maps onto the
//! existing [`MSSpectrum::retain_peaks`] / [`MSChromatogram::retain_peaks`].
//!
//! The source templates over any peak container with positions, intensities and
//! the three data-array lists. The [`PeakContainer`] trait carries exactly that
//! surface for [`MSSpectrum`] and [`MSChromatogram`], so the helpers are generic
//! where the source is generic. Operations whose cost scales with the peak count
//! are checked against [`SpectrumHelperLimits`] before anything is mutated; a
//! failure leaves the container unchanged.

use super::{ChromatogramPeak, DataArray, MSChromatogram, MSSpectrum, Peak1D};
use crate::error::{Error, Result};
use std::cmp::Ordering;

/// Peak container surface shared by [`MSSpectrum`] and [`MSChromatogram`].
///
/// This is the trait counterpart of the source's `PeakContainerT` template
/// parameter: one coordinate per peak (`m/z` or retention time), an `f32`
/// intensity, and the float/integer/string data-array lists. Implementors are
/// the two kernel containers; the trait is public so callers can write their own
/// generic helpers over both.
pub trait PeakContainer: Default {
    /// Peak value type (`Peak1D` or `ChromatogramPeak`).
    type Peak: Copy;

    /// Borrowed peaks in storage order.
    fn peaks(&self) -> &[Self::Peak];

    /// Mutable peak storage; replacing it does not touch data arrays.
    fn peaks_mut(&mut self) -> &mut Vec<Self::Peak>;

    /// Coordinate of a peak: m/z for spectra, retention time in seconds for
    /// chromatograms (source `getPos`).
    fn position(peak: &Self::Peak) -> f64;

    /// Intensity of a peak (source `getIntensity`).
    fn intensity(peak: &Self::Peak) -> f32;

    /// Overwrite a peak's intensity (source `setIntensity`).
    fn set_intensity(peak: &mut Self::Peak, intensity: f32);

    /// Construct a peak from coordinate and intensity (source
    /// `PeakType(position, intensity)`).
    fn new_peak(position: f64, intensity: f32) -> Self::Peak;

    /// True when any float, integer or string data array is attached, even one
    /// with no entries. This mirrors the source test on the array *lists*.
    fn has_data_arrays(&self) -> bool;

    /// Remove every data array, leaving peaks and metadata alone.
    fn clear_data_arrays(&mut self);
}

macro_rules! peak_container_impl {
    ($container:ty, $peak:ty, $position:ident) => {
        impl PeakContainer for $container {
            type Peak = $peak;

            fn peaks(&self) -> &[Self::Peak] {
                &self.peaks
            }

            fn peaks_mut(&mut self) -> &mut Vec<Self::Peak> {
                &mut self.peaks
            }

            fn position(peak: &Self::Peak) -> f64 {
                peak.$position
            }

            fn intensity(peak: &Self::Peak) -> f32 {
                peak.intensity
            }

            fn set_intensity(peak: &mut Self::Peak, intensity: f32) {
                peak.intensity = intensity;
            }

            fn new_peak(position: f64, intensity: f32) -> Self::Peak {
                <$peak>::new(position, intensity)
            }

            fn has_data_arrays(&self) -> bool {
                !self.float_data_arrays.is_empty()
                    || !self.integer_data_arrays.is_empty()
                    || !self.string_data_arrays.is_empty()
            }

            fn clear_data_arrays(&mut self) {
                self.float_data_arrays.clear();
                self.integer_data_arrays.clear();
                self.string_data_arrays.clear();
            }
        }
    };
}

peak_container_impl!(MSSpectrum, Peak1D, mz);
peak_container_impl!(MSChromatogram, ChromatogramPeak, rt);

/// Per-call ceilings for the helpers that scale with the peak count.
///
/// `max_peaks` bounds the number of peaks accepted; `max_work` bounds the
/// logical visit count charged before the operation starts. Sorting charges a
/// conservative `8 * n * bit_length(n)` allowance covering the position sort and
/// the per-group median sorts, plus `2 * n` for the merge pass; intensity
/// rebasing charges `3 * n`. The defaults admit the ten-million-peak ceiling
/// under the work ceiling. These are logical bounds, not allocator measurements.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpectrumHelperLimits {
    /// Maximum number of peaks accepted by one call.
    pub max_peaks: usize,
    /// Maximum logical work units charged by one call.
    pub max_work: usize,
}

impl Default for SpectrumHelperLimits {
    fn default() -> Self {
        Self {
            max_peaks: 10_000_000,
            max_work: 4_000_000_000,
        }
    }
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

fn limit() -> Error {
    invalid("spectrum helper resource limit exceeded")
}

fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}

fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}

fn bit_length(value: usize) -> usize {
    (usize::BITS - value.leading_zeros()) as usize
}

/// Check the peak-count and work ceilings before any allocation or mutation.
fn preflight(len: usize, limits: SpectrumHelperLimits, sorting: bool) -> Result<()> {
    if len > limits.max_peaks {
        return Err(limit());
    }
    let work = if sorting {
        add(mul(2, len)?, mul(mul(8, len)?, bit_length(len))?)?
    } else {
        mul(3, len)?
    };
    if work > limits.max_work {
        return Err(limit());
    }
    Ok(())
}

/// Index of the first data array with the given name.
///
/// Source `getDataArrayByName` returns an iterator, or the end iterator when no
/// array carries `name`; this returns `Some(index)` or `None`. One generic
/// function covers the source's mutable and const overloads for float, integer
/// and string arrays of spectra and chromatograms. Names compare exactly. The
/// cost is one linear scan over the array list, which the caller already owns.
pub fn data_array_index_by_name<T>(arrays: &[DataArray<T>], name: &str) -> Option<usize> {
    arrays.iter().position(|array| array.name == name)
}

/// Borrow the first data array with the given name, or `None`.
///
/// Counterpart of the source const `getDataArrayByName` overload; see
/// [`data_array_index_by_name`] for the index form.
pub fn data_array_by_name<'a, T>(
    arrays: &'a [DataArray<T>],
    name: &str,
) -> Option<&'a DataArray<T>> {
    data_array_index_by_name(arrays, name).map(|index| &arrays[index])
}

/// Mutably borrow the first data array with the given name, or `None`.
///
/// Counterpart of the source mutable `getDataArrayByName` overload; see
/// [`data_array_index_by_name`] for the index form.
pub fn data_array_by_name_mut<'a, T>(
    arrays: &'a mut [DataArray<T>],
    name: &str,
) -> Option<&'a mut DataArray<T>> {
    data_array_index_by_name(arrays, name).map(move |index| &mut arrays[index])
}

/// Shift every intensity so the minimum becomes zero, using default limits.
///
/// See [`subtract_minimum_intensity_with_limits`].
///
/// # Errors
///
/// As for [`subtract_minimum_intensity_with_limits`] with
/// [`SpectrumHelperLimits::default`].
pub fn subtract_minimum_intensity<C: PeakContainer>(container: &mut C) -> Result<()> {
    subtract_minimum_intensity_with_limits(container, SpectrumHelperLimits::default())
}

/// Shift every intensity so the minimum becomes zero.
///
/// The source finds the first minimum intensity, forms `rebase = -minimum` as a
/// `double`, and stores `float(intensity + rebase)` for every peak; the same
/// `f32 -> f64 -> f32` arithmetic is used here, so a negative minimum raises all
/// intensities and a positive minimum lowers them. An empty container is left
/// unchanged, as in the source.
///
/// The source notes that data arrays are not updated. Intensities are not stored
/// in data arrays, so nothing is lost; the arrays are left untouched here too.
/// Peak order and coordinates are unchanged.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `container` exceeds `limits`, when any
/// intensity is not finite, or when a rebased intensity would not be
/// representable as a finite `f32`. The source performs none of these checks;
/// all of them run before the first intensity is written, so a failure leaves
/// the container unchanged.
pub fn subtract_minimum_intensity_with_limits<C: PeakContainer>(
    container: &mut C,
    limits: SpectrumHelperLimits,
) -> Result<()> {
    let peaks = container.peaks();
    preflight(peaks.len(), limits, false)?;
    let Some(first) = peaks.first() else {
        return Ok(());
    };
    let mut minimum = C::intensity(first);
    for peak in peaks {
        let intensity = C::intensity(peak);
        if !intensity.is_finite() {
            return Err(invalid("peak intensity must be finite"));
        }
        if intensity < minimum {
            minimum = intensity;
        }
    }
    let rebase = -f64::from(minimum);
    if peaks
        .iter()
        .any(|peak| !rebased(C::intensity(peak), rebase).is_finite())
    {
        return Err(invalid("rebased peak intensity is not a finite f32"));
    }
    for peak in container.peaks_mut().iter_mut() {
        let intensity = rebased(C::intensity(peak), rebase);
        C::set_intensity(peak, intensity);
    }
    Ok(())
}

/// Source expression `float(intensity + rebase)` with `rebase` a `double`.
fn rebased(intensity: f32, rebase: f64) -> f32 {
    (f64::from(intensity) + rebase) as f32
}

/// How intensities of peaks sharing a position are combined.
///
/// Source `IntensityAveragingMethod`, declared in this order as `MEDIAN, MEAN,
/// SUM, MIN, MAX`. The source default argument is `MEDIAN`, which is the
/// `Default` here. See [`make_peak_position_unique`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum IntensityAveragingMethod {
    /// Sorted median; an even count averages the two middle values.
    #[default]
    Median,
    /// Arithmetic mean of the group.
    Mean,
    /// Left-to-right sum of the group.
    Sum,
    /// Smallest intensity in the group.
    Min,
    /// Largest intensity in the group.
    Max,
}

/// Options for [`make_peak_position_unique_with`].
///
/// The native default keeps the container's metadata and refuses to run while
/// data arrays are attached. [`UniquePositionOptions::source`] reproduces the
/// source, which discards both.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UniquePositionOptions {
    /// Drop attached float/integer/string data arrays instead of refusing.
    /// The source drops them after logging a warning; per-peak annotations
    /// cannot be merged, so the native default is to refuse.
    pub discard_data_arrays: bool,
    /// Replace the whole container with a default-constructed one holding only
    /// the merged peaks. The source does this through `std::swap` with a fresh
    /// container, losing RT, MS level, name, settings and everything else; the
    /// native default keeps all non-peak fields.
    pub reset_metadata: bool,
    /// Per-call ceilings checked before any mutation.
    pub limits: SpectrumHelperLimits,
}

impl UniquePositionOptions {
    /// Source behaviour: data arrays are discarded and every other field is
    /// reset to its default, exactly as the source's swap with a fresh container.
    pub fn source() -> Self {
        Self {
            discard_data_arrays: true,
            reset_metadata: true,
            limits: SpectrumHelperLimits::default(),
        }
    }
}

/// Merge peaks sharing a position, keeping metadata and default limits.
///
/// Equivalent to [`make_peak_position_unique_with`] using
/// [`UniquePositionOptions::default`]: metadata is kept, attached data arrays
/// are refused, and the default [`SpectrumHelperLimits`] apply.
///
/// # Errors
///
/// As for [`make_peak_position_unique_with`].
///
/// ```rust
/// use openms::kernel::spectrum_helper::{make_peak_position_unique, IntensityAveragingMethod};
/// use openms::kernel::{MSSpectrum, Peak1D};
///
/// let mut spectrum = MSSpectrum::from_peaks(vec![
///     Peak1D::new(1.0, 1.0),
///     Peak1D::new(2.0, 4.0),
///     Peak1D::new(2.0, 8.0),
///     Peak1D::new(2.0, 10.0),
/// ]);
/// spectrum.rt = 12.5;
/// make_peak_position_unique(&mut spectrum, IntensityAveragingMethod::Sum).unwrap();
/// assert_eq!(spectrum.peaks, vec![Peak1D::new(1.0, 1.0), Peak1D::new(2.0, 22.0)]);
/// assert_eq!(spectrum.rt, 12.5); // native: metadata survives (the source resets it)
/// ```
pub fn make_peak_position_unique<C: PeakContainer>(
    container: &mut C,
    method: IntensityAveragingMethod,
) -> Result<()> {
    make_peak_position_unique_with(container, method, UniquePositionOptions::default())
}

/// Make peak positions unique.
///
/// A peak container may contain multiple peaks with the same position, i.e.
/// either a spectrum containing peaks with the same m/z position, or a
/// chromatogram containing peaks with identical RT position. One scenario where
/// this happens is when multiple spectra are merged into a single one. The
/// method combines peaks with the same position into a single one whose
/// intensity is determined by `method`.
///
/// The source algorithm is reproduced exactly: the peaks are stably sorted by
/// position, walked once, and a new group starts whenever a position is
/// strictly greater than the current one (so `-0.0` and `0.0` share a group).
/// Each group's intensities are widened to `f64` in storage order; the median
/// sorts them and averages the two middle values for an even count, the mean is
/// the sum divided by the count, the sum is the left-to-right `f64`
/// accumulation from zero, and min/max are the extreme values. The `f64` result
/// is narrowed to the `f32` peak intensity, as the source's peak constructor
/// does. Merged peaks are emitted in ascending position order. An empty
/// container is returned unchanged before any policy check, as in the source.
///
/// # Source behaviour and native differences
///
/// The source warns that data arrays are ignored, then finishes with
/// `std::swap(p_new, p)` against a *default-constructed* container. The result
/// therefore loses not only the data arrays but every other field: RT, MS level,
/// name, native ID, precursors, settings and metadata. With the default
/// [`UniquePositionOptions`] this port keeps all non-peak fields and returns an
/// error when any data array is attached; `discard_data_arrays` drops the
/// arrays, and `reset_metadata` reproduces the source's reset.
/// [`UniquePositionOptions::source`] selects both. No warning is logged; the
/// refusal or the explicit option replaces it.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the peak count or the sorting work
/// exceeds `options.limits`, when data arrays are attached and
/// `options.discard_data_arrays` is false, when any position or intensity is not
/// finite (the source would sort with an unordered `NaN`), or when a merged
/// intensity is not a finite `f32` (for example a `Sum` overflowing `f32`). The
/// merged peaks are built in a temporary and committed only at the end, so any
/// error leaves the container unchanged.
pub fn make_peak_position_unique_with<C: PeakContainer>(
    container: &mut C,
    method: IntensityAveragingMethod,
    options: UniquePositionOptions,
) -> Result<()> {
    let peaks = container.peaks();
    preflight(peaks.len(), options.limits, true)?;
    if peaks.is_empty() {
        return Ok(());
    }
    if container.has_data_arrays() && !options.discard_data_arrays {
        return Err(invalid(
            "data arrays are attached; makePeakPositionUnique cannot merge them \
             (set UniquePositionOptions::discard_data_arrays to drop them as the source does)",
        ));
    }
    for peak in peaks {
        if !C::position(peak).is_finite() {
            return Err(invalid("peak position must be finite"));
        }
        if !C::intensity(peak).is_finite() {
            return Err(invalid("peak intensity must be finite"));
        }
    }

    let mut sorted = Vec::new();
    sorted.try_reserve_exact(peaks.len()).map_err(|_| limit())?;
    sorted.extend_from_slice(peaks);
    // Stable, like the source std::stable_sort, so intensities keep storage
    // order within a group; the sum's accumulation order depends on it.
    sorted.sort_by(|a, b| {
        C::position(a)
            .partial_cmp(&C::position(b))
            .unwrap_or(Ordering::Equal)
    });

    let mut merged = Vec::new();
    let mut group = Vec::new();
    let mut current = C::position(&sorted[0]);
    for peak in &sorted {
        let position = C::position(peak);
        if position > current {
            merged.push(C::new_peak(current, combine(method, &mut group)?));
            current = position;
            group.clear();
        }
        group.push(f64::from(C::intensity(peak)));
    }
    merged.push(C::new_peak(current, combine(method, &mut group)?));

    if options.reset_metadata {
        let mut fresh = C::default();
        *fresh.peaks_mut() = merged;
        *container = fresh;
    } else {
        *container.peaks_mut() = merged;
        if options.discard_data_arrays {
            container.clear_data_arrays();
        }
    }
    Ok(())
}

/// Combine one non-empty group with the source `Math::median/mean/sum` and
/// `std::min_element/max_element` semantics, then narrow to `f32`.
fn combine(method: IntensityAveragingMethod, values: &mut [f64]) -> Result<f32> {
    let count = values.len();
    let sum = || values.iter().fold(0.0_f64, |sum, value| sum + value);
    let combined = match method {
        IntensityAveragingMethod::Median => {
            values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
            if count % 2 == 0 {
                (values[count / 2 - 1] + values[count / 2]) / 2.0
            } else {
                values[(count - 1) / 2]
            }
        }
        IntensityAveragingMethod::Mean => sum() / count as f64,
        IntensityAveragingMethod::Sum => sum(),
        IntensityAveragingMethod::Min => values.iter().copied().fold(f64::INFINITY, f64::min),
        IntensityAveragingMethod::Max => values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    };
    let narrowed = combined as f32;
    if narrowed.is_finite() {
        Ok(narrowed)
    } else {
        Err(invalid("merged peak intensity is not a finite f32"))
    }
}

/// Copy only the metadata of `input` onto `output`; peak data is not copied.
///
/// Source `copySpectrumMeta`: when `clear_spectrum` is true the output is first
/// cleared with `clear(true)` (all peaks, data arrays and metadata), then the
/// spectrum settings, RT, drift time and unit, MS level and name are assigned
/// from the input. When it is false the output keeps its peaks and data arrays
/// and only the metadata is overwritten. Here every field of [`MSSpectrum`] other
/// than `peaks` and the three data-array lists is cloned from `input`; with
/// `clear_spectrum` those four are emptied, otherwise they are retained.
///
/// The native `peptide_identifications` field, which the source keeps outside
/// `MSSpectrum`, is treated as metadata and copied. The source's explicit
/// `setDriftTime`/`setDriftTimeUnit` calls are covered by the same clone, since
/// [`MSSpectrum::drift_time`] and [`MSSpectrum::drift_time_unit`] are ordinary
/// fields here. No ceiling applies: the call clones metadata only, which is
/// outside the checked-operation budgets like ordinary `Clone`.
pub fn copy_spectrum_meta(input: &MSSpectrum, output: &mut MSSpectrum, clear_spectrum: bool) {
    let (peaks, float_data_arrays, integer_data_arrays, string_data_arrays) = if clear_spectrum {
        (Vec::new(), Vec::new(), Vec::new(), Vec::new())
    } else {
        (
            std::mem::take(&mut output.peaks),
            std::mem::take(&mut output.float_data_arrays),
            std::mem::take(&mut output.integer_data_arrays),
            std::mem::take(&mut output.string_data_arrays),
        )
    };
    *output = MSSpectrum {
        peaks,
        float_data_arrays,
        integer_data_arrays,
        string_data_arrays,
        ..input.clone()
    };
}
