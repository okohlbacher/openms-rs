// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A binned counter over a closed value range.
//!
//! Port of `src/openms/include/OpenMS/MATH/STATISTICS/Histogram.h`; the
//! accompanying `Histogram.cpp` is an empty translation unit. See
//! `docs/HISTOGRAM_SUPPORT.md`.
//!
//! The source is `template <typename ValueType = UInt, typename BinSizeType =
//! double> class Histogram`, so the bin contents and the coordinate type are
//! independent. This port fixes both to `f64`. That keeps every instantiation
//! the library uses expressible — including the `Histogram<float, float>` of
//! the class test, whose bins hold fractional increments — and removes the
//! integer truncation that
//! [`crate::math::histogram::Histogram::apply_log_transformation`] suffers
//! under the default `ValueType = UInt`.
//!
//! The binning convention is what the boundary cases turn on. Bins are
//! half-open, `[left, right)`, and are laid out from `min` with width
//! `bin_size`; the number of bins is `ceil((max - min) / bin_size)`, so the
//! last bin usually extends past `max`. `max` itself is **in range** and always
//! lands in the last bin, whatever the arithmetic would otherwise say; that
//! special case is why
//! [`crate::math::histogram::Histogram::right_border_of_bin`] reports the
//! next representable value above `max` for the last bin rather than a
//! multiple of `bin_size`.

use crate::{Error, Result};

/// Maximum number of bins one histogram may hold.
///
/// `ceil((max - min) / bin_size)` is entirely caller-controlled, and the source
/// hands the result straight to `std::vector::resize`. Native-only guard.
pub const MAX_BINS: usize = 16_777_216;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

fn out_of_range(message: &str) -> Error {
    Error::InvalidRange(message.to_string())
}

/// Refuse a non-finite bin increment before any bin is touched.
fn check_increment(increment: f64) -> Result<()> {
    if !increment.is_finite() {
        return Err(bad("histogram increment must be finite"));
    }
    Ok(())
}

/// The next representable `f64` above `value`, towards `+inf`.
///
/// The source calls `std::nextafter(maxBound(), maxBound() + 1)`. `f64::next_up`
/// is not available at this crate's minimum supported Rust version, so the step
/// is taken on the bit pattern, which is the same operation for every finite
/// input.
fn next_above(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        // Covers -0.0 as well; the next value above either zero is the
        // smallest positive subnormal.
        return f64::from_bits(1);
    }
    if value > 0.0 {
        f64::from_bits(value.to_bits() + 1)
    } else {
        f64::from_bits(value.to_bits() - 1)
    }
}

/// A histogram over the closed range `[min, max]` with fixed-width bins.
///
/// Equality compares the bounds, the bin width and every bin, as the source's
/// `operator==` does.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Histogram {
    min: f64,
    max: f64,
    bin_size: f64,
    bins: Vec<f64>,
}

impl Histogram {
    /// An empty histogram with zero bounds, zero bin width and no bins.
    ///
    /// This is the source's default constructor, kept because the class test
    /// constructs one and because it is the natural target of an assignment.
    /// It has no bins, so every accessor on it reports an error rather than
    /// dereferencing past the end of an empty vector as the source does.
    pub fn new() -> Self {
        Self::default()
    }

    /// A histogram over `[min, max]` with bins of width `bin_size`.
    ///
    /// The bin count is `ceil((max - min) / bin_size)`, except that `max ==
    /// min` yields exactly one bin.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `bin_size` is not positive, which
    /// is the source's `Exception::OutOfRange`, and additionally when any
    /// argument is not finite, when `max < min`, or when the resulting bin
    /// count exceeds [`MAX_BINS`]. The source checks only the bin width: an
    /// inverted range makes it take the ceiling of a negative quotient and
    /// convert it to an unsigned size, and a narrow bin width lets the bin
    /// count grow until the allocation fails.
    pub fn with_bounds(min: f64, max: f64, bin_size: f64) -> Result<Self> {
        let bins = vec![0.0; Self::bin_count(min, max, bin_size)?];
        Ok(Self {
            min,
            max,
            bin_size,
            bins,
        })
    }

    /// A histogram over `[min, max]`, filled from `values`.
    ///
    /// Equivalent to [`Histogram::with_bounds`] followed by one
    /// [`Histogram::inc`] per value, in slice order.
    ///
    /// # Errors
    ///
    /// As [`Histogram::with_bounds`], plus [`Error::InvalidRange`] when any
    /// value lies outside `[min, max]` or is not finite. The source propagates
    /// the same `Exception::OutOfRange` out of the constructor, leaving the
    /// half-filled object unreachable.
    pub fn from_values(values: &[f64], min: f64, max: f64, bin_size: f64) -> Result<Self> {
        let mut histogram = Self::with_bounds(min, max, bin_size)?;
        for &value in values {
            histogram.inc(value)?;
        }
        Ok(histogram)
    }

    /// Number of bins spanning `[min, max]` at `bin_size`, with every guard.
    fn bin_count(min: f64, max: f64, bin_size: f64) -> Result<usize> {
        if !min.is_finite() || !max.is_finite() || !bin_size.is_finite() {
            return Err(bad("histogram bounds and bin size must be finite"));
        }
        if bin_size <= 0.0 {
            return Err(bad("histogram bin size must be positive"));
        }
        if max < min {
            return Err(bad(
                "histogram upper bound must not be below the lower bound",
            ));
        }
        if max == min {
            // The source's explicit special case: one bin, not zero.
            return Ok(1);
        }
        let count = ((max - min) / bin_size).ceil();
        if !count.is_finite() || count > MAX_BINS as f64 {
            return Err(bad("histogram bin count exceeds the supported maximum"));
        }
        // A quotient that underflows to zero would leave the source with an
        // empty bin vector and a range that nothing can be counted into.
        Ok((count as usize).max(1))
    }

    /// The lower bound, inclusive.
    pub fn min_bound(&self) -> f64 {
        self.min
    }

    /// The upper bound, inclusive.
    pub fn max_bound(&self) -> f64 {
        self.max
    }

    /// The bin width.
    pub fn bin_size(&self) -> f64 {
        self.bin_size
    }

    /// The number of bins.
    pub fn len(&self) -> usize {
        self.bins.len()
    }

    /// Whether the histogram has no bins, which only a default-constructed one
    /// does.
    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    /// The highest bin count, or `None` when there are no bins.
    ///
    /// The source returns `*std::max_element(...)` unconditionally, which
    /// dereferences the end iterator of a default-constructed histogram.
    pub fn max_value(&self) -> Option<f64> {
        self.bins.iter().copied().reduce(f64::max)
    }

    /// The lowest bin count, or `None` when there are no bins.
    ///
    /// See [`Histogram::max_value`] for the source's behaviour without bins.
    pub fn min_value(&self) -> Option<f64> {
        self.bins.iter().copied().reduce(f64::min)
    }

    /// The count in bin `index`.
    ///
    /// The source spells this `operator[]`, which throws on an invalid index.
    /// A Rust `Index` implementation would have to panic instead, so this is a
    /// named method returning an error.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for an index at or beyond
    /// [`Histogram::len`], matching the source's `Exception::IndexOverflow`.
    pub fn bin(&self, index: usize) -> Result<f64> {
        self.bins
            .get(index)
            .copied()
            .ok_or_else(|| out_of_range("histogram bin index out of range"))
    }

    /// All bins, in order.
    pub fn bins(&self) -> &[f64] {
        &self.bins
    }

    /// Iterator over the bins, in order.
    ///
    /// The source's `begin()`/`end()` pair over its `ConstIterator`.
    pub fn iter(&self) -> std::slice::Iter<'_, f64> {
        self.bins.iter()
    }

    /// The centre coordinate of bin `index`, `min + (index + 0.5) * bin_size`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for an index at or beyond
    /// [`Histogram::len`].
    pub fn center_of_bin(&self, index: usize) -> Result<f64> {
        self.check_index(index)?;
        Ok(self.min + (index as f64 + 0.5) * self.bin_size)
    }

    /// The leftmost coordinate belonging to bin `index`, `min + index *
    /// bin_size`. The interval is closed on this side.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for an index at or beyond
    /// [`Histogram::len`].
    pub fn left_border_of_bin(&self, index: usize) -> Result<f64> {
        self.check_index(index)?;
        Ok(self.min + index as f64 * self.bin_size)
    }

    /// The first coordinate to the right of bin `index` that is no longer part
    /// of it. The interval is open on this side.
    ///
    /// For the last bin this is the next representable value above
    /// [`Histogram::max_bound`], not `min + (index + 1) * bin_size`: the last
    /// bin is special because [`Histogram::value_to_bin`] routes `max` into it
    /// regardless of the bin arithmetic.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] for an index at or beyond
    /// [`Histogram::len`].
    pub fn right_border_of_bin(&self, index: usize) -> Result<f64> {
        self.check_index(index)?;
        if index + 1 == self.bins.len() {
            return Ok(next_above(self.max));
        }
        Ok(self.min + (index as f64 + 1.0) * self.bin_size)
    }

    /// The count in the bin that `value` falls into.
    ///
    /// # Errors
    ///
    /// As [`Histogram::value_to_bin`].
    pub fn bin_value(&self, value: f64) -> Result<f64> {
        let index = self.value_to_bin(value)?;
        self.bin(index)
    }

    /// Add one to the bin that `value` falls into, returning that bin's index.
    ///
    /// # Errors
    ///
    /// As [`Histogram::value_to_bin`].
    pub fn inc(&mut self, value: f64) -> Result<usize> {
        self.inc_by(value, 1.0)
    }

    /// Add `increment` to the bin that `value` falls into, returning that
    /// bin's index.
    ///
    /// The source's `inc(val, increment = 1)`. The increment is not required to
    /// be positive or integral; the class test uses fractional counts.
    ///
    /// # Errors
    ///
    /// As [`Histogram::value_to_bin`], plus [`Error::InvalidValue`] for a
    /// non-finite increment, which the source would store and which would then
    /// poison [`Histogram::min_value`] and [`Histogram::max_value`]. The bin is
    /// left untouched when either check fails.
    pub fn inc_by(&mut self, value: f64, increment: f64) -> Result<usize> {
        check_increment(increment)?;
        let index = self.value_to_bin(value)?;
        self.bins[index] += increment;
        Ok(index)
    }

    /// Add `increment` to every bin below the bin `value` falls into,
    /// returning that bin's index.
    ///
    /// With `inclusive`, the bin of `value` is incremented too.
    ///
    /// # Errors
    ///
    /// As [`Histogram::inc_by`]. No bin is touched when either check fails.
    pub fn inc_until(&mut self, value: f64, inclusive: bool, increment: f64) -> Result<usize> {
        check_increment(increment)?;
        let index = self.value_to_bin(value)?;
        for bin in self.bins.iter_mut().take(index) {
            *bin += increment;
        }
        if inclusive {
            self.bins[index] += increment;
        }
        Ok(index)
    }

    /// Add `increment` to every bin above the bin `value` falls into,
    /// returning that bin's index.
    ///
    /// With `inclusive`, the bin of `value` is incremented too.
    ///
    /// # Errors
    ///
    /// As [`Histogram::inc_by`]. No bin is touched when either check fails.
    pub fn inc_from(&mut self, value: f64, inclusive: bool, increment: f64) -> Result<usize> {
        check_increment(increment)?;
        let index = self.value_to_bin(value)?;
        for bin in self.bins.iter_mut().skip(index + 1) {
            *bin += increment;
        }
        if inclusive {
            self.bins[index] += increment;
        }
        Ok(index)
    }

    /// Accumulate a cumulative histogram from `values`.
    ///
    /// The source's static `getCumulativeHistogram`, which takes the histogram
    /// by reference and adds to whatever it already holds. With `complement`,
    /// each value raises every bin below it through
    /// [`Histogram::inc_until`]; otherwise it raises every bin above it through
    /// [`Histogram::inc_from`]. `inclusive` decides whether the value's own bin
    /// is raised.
    ///
    /// # Errors
    ///
    /// As [`Histogram::value_to_bin`], on the first offending value. Values
    /// already applied stay applied, exactly as in the source, which has no
    /// rollback either; check the range first when that matters.
    pub fn add_cumulative(
        &mut self,
        values: &[f64],
        complement: bool,
        inclusive: bool,
    ) -> Result<()> {
        for &value in values {
            if complement {
                self.inc_until(value, inclusive, 1.0)?;
            } else {
                self.inc_from(value, inclusive, 1.0)?;
            }
        }
        Ok(())
    }

    /// Replace the range and bin width, discarding all counts.
    ///
    /// # Errors
    ///
    /// As [`Histogram::with_bounds`]. The new bins are built before anything is
    /// replaced, so a refused reset leaves the histogram exactly as it was; the
    /// source clears its bins first and then throws, leaving an object whose
    /// bounds have changed and whose bins are gone.
    pub fn reset(&mut self, min: f64, max: f64, bin_size: f64) -> Result<()> {
        let bins = vec![0.0; Self::bin_count(min, max, bin_size)?];
        self.min = min;
        self.max = max;
        self.bin_size = bin_size;
        self.bins = bins;
        Ok(())
    }

    /// Transform every bin by `f(x) = multiplier * ln(x + 1)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when any bin is at or below `-1`, where
    /// the logarithm is undefined; the source evaluates it anyway and stores
    /// NaN. Bins are transformed into a temporary and committed together, so a
    /// refusal leaves every bin unchanged.
    ///
    /// The source casts each result back to its `ValueType`, so its default
    /// `UInt` instantiation truncates the transformed value to an integer. This
    /// port stores `f64` and does not truncate; the class test instantiates
    /// `Histogram<float, float>` and therefore also sees untruncated values.
    pub fn apply_log_transformation(&mut self, multiplier: f64) -> Result<()> {
        let mut transformed = Vec::with_capacity(self.bins.len());
        for &bin in &self.bins {
            let argument = bin + 1.0;
            if argument <= 0.0 || !argument.is_finite() {
                return Err(bad("log transformation needs every bin above -1"));
            }
            transformed.push(multiplier * argument.ln());
        }
        self.bins = transformed;
        Ok(())
    }

    /// The index of the bin `value` belongs to.
    ///
    /// `floor((value - min) / bin_size)`, except that `value == max` is routed
    /// into the last bin: the range is closed at both ends, but the bins are
    /// half-open, so without that case the upper bound would have no bin.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidRange`] when `value` is outside `[min, max]`,
    /// matching the source's `Exception::OutOfRange`, and additionally when
    /// `value` is NaN or the histogram has no bins. The source's two bound
    /// comparisons are both false for a NaN, which then reaches an unsigned
    /// conversion of `floor(NaN)`.
    pub fn value_to_bin(&self, value: f64) -> Result<usize> {
        if self.bins.is_empty() {
            return Err(out_of_range("histogram has no bins"));
        }
        if value.is_nan() {
            return Err(out_of_range("histogram value must not be NaN"));
        }
        if value < self.min || value > self.max {
            return Err(out_of_range("value outside the histogram range"));
        }
        if value == self.max {
            return Ok(self.bins.len() - 1);
        }
        let position = ((value - self.min) / self.bin_size).floor();
        if position < 0.0 {
            return Err(out_of_range("value outside the histogram range"));
        }
        // `position` is below the bin count for every value below `max`; the
        // clamp only catches a division that rounded up onto the boundary.
        Ok((position as usize).min(self.bins.len() - 1))
    }

    fn check_index(&self, index: usize) -> Result<()> {
        if index >= self.bins.len() {
            return Err(out_of_range("histogram bin index out of range"));
        }
        Ok(())
    }
}

impl<'a> IntoIterator for &'a Histogram {
    type Item = &'a f64;
    type IntoIter = std::slice::Iter<'a, f64>;

    fn into_iter(self) -> Self::IntoIter {
        self.bins.iter()
    }
}

impl std::fmt::Display for Histogram {
    /// The source's `operator<<`: one line per bin, its centre and its count
    /// separated by a tab.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (index, bin) in self.bins.iter().enumerate() {
            let center = self.min + (index as f64 + 0.5) * self.bin_size;
            writeln!(f, "{center}\t{bin}")?;
        }
        Ok(())
    }
}
