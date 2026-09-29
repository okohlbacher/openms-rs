// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Linearly interpolated values and derivatives from equally spaced samples.
//!
//! Port of `src/openms/include/OpenMS/ML/INTERPOLATION/LinearInterpolation.h`,
//! a header-only template (`template <typename Key = double, typename Value =
//! Key>`) with no translation unit. The library instantiates it only as
//! `<double, double>` (through `InterpolationModel`), so
//! [`LinearInterpolation`](crate::math::ml::linear_interpolation::LinearInterpolation) is not generic: key and value are `f64`. The class
//! test also runs a `<float, double>` instantiation; its sections are ported
//! with `f64` keys, which changes nothing for its dyadic sample positions. See
//! `docs/FEATUREFINDER_MODELS_SUPPORT.md`.
//!
//! This module depends on nothing above `math`, as the `ml` root requires.

/// Linearly interpolated values (and derivatives) from discrete samples.
///
/// The samples ("data") are the y-coordinates at the "inside" positions
/// `0, 1, ..., data.len() - 1`. Values beyond that range are implicitly zero.
/// The x-axis can be scaled and shifted by an affine mapping between "outside"
/// and "inside" coordinates, given either by [`LinearInterpolation::set_scale`]
/// and [`LinearInterpolation::set_offset`] or by
/// [`LinearInterpolation::set_mapping`]. The default is the identity mapping
/// (scale 1, offset 0).
///
/// [`LinearInterpolation::value`] and [`LinearInterpolation::derivative`] sample
/// the interpolated data and its derivative at an outside position.
#[derive(Clone, Debug, PartialEq)]
pub struct LinearInterpolation {
    scale: f64,
    offset: f64,
    inside: f64,
    outside: f64,
    data: Vec<f64>,
}

impl Default for LinearInterpolation {
    /// Scale 1, offset 0, no data: the source's default arguments.
    fn default() -> Self {
        Self::new(1.0, 0.0)
    }
}

/// Split an inside position into its integer part as an index and its
/// fractional part, as the source's `std::modf` and `DiffType(left_key)` do.
///
/// `None` stands for an integer part that no index can hold (an infinity or a
/// magnitude beyond `isize`), where the source converts out of range, which is
/// undefined behaviour; every caller treats it as "far outside the data".
fn split(pos: f64) -> Option<(isize, f64)> {
    let whole = pos.trunc();
    // isize::MAX as f64 rounds up to 2^63, so the strict bound excludes it.
    if !(whole > isize::MIN as f64 && whole < isize::MAX as f64) {
        return None;
    }
    Some((whole as isize, pos - whole))
}

impl LinearInterpolation {
    /// An interpolation with the given `scale` and `offset` and no data.
    ///
    /// `scale` is applied to the arguments of [`Self::value`] and
    /// [`Self::derivative`] before the lookup in the data; `offset` is
    /// subtracted before everything else. The inside and outside reference
    /// points start at zero.
    pub fn new(scale: f64, offset: f64) -> Self {
        Self {
            scale,
            offset,
            inside: 0.0,
            outside: 0.0,
            data: Vec::new(),
        }
    }

    /// The interpolated value at outside position `arg_pos`.
    ///
    /// Between two samples the value is their linear interpolation; within one
    /// inside unit beyond either end it falls linearly to zero, and further out
    /// it is zero.
    ///
    /// Native detail: a NaN position yields NaN; an infinite position, or one
    /// whose inside index exceeds `isize`, yields zero; and empty data yields
    /// zero everywhere. The source reads `data[0]` of an empty container for
    /// positions in `(-1, 0)` and converts a non-finite index to an integer,
    /// both undefined behaviour.
    pub fn value(&self, arg_pos: f64) -> f64 {
        let pos = self.key2index(arg_pos);
        if pos.is_nan() {
            return f64::NAN;
        }
        let Some((left, frac)) = split(pos) else {
            return 0.0;
        };
        if self.data.is_empty() {
            return 0.0;
        }
        if pos < 0.0 {
            if left != 0 {
                return 0.0;
            }
            return self.data[0] * (1.0 + frac);
        }
        // `left` is non-negative here and data is non-empty.
        let back = self.data.len() - 1;
        let left = left as usize;
        if left >= back {
            if left != back {
                return 0.0;
            }
            return self.data[left] * (1.0 - frac);
        }
        self.data[left + 1] * frac + self.data[left] * (1.0 - frac)
    }

    /// Linear resampling: split `arg_value` between the two samples around
    /// outside position `arg_pos` and add the parts to them.
    ///
    /// The weights are those [`Self::value`] would use at `arg_pos`, so a value
    /// added one inside unit or less beyond either end lands partly on the end
    /// sample, and a value added further out is dropped.
    ///
    /// Native detail: a NaN, infinite or out-of-`isize` position, and any
    /// position on empty data, adds nothing; the source is undefined there, as
    /// for [`Self::value`].
    pub fn add_value(&mut self, arg_pos: f64, arg_value: f64) {
        let pos = self.key2index(arg_pos);
        if pos.is_nan() || self.data.is_empty() {
            return;
        }
        let Some((left, frac)) = split(pos) else {
            return;
        };
        if pos < 0.0 {
            if left == 0 {
                self.data[0] += (1.0 + frac) * arg_value;
            }
            return;
        }
        let back = self.data.len() - 1;
        let left = left as usize;
        if left >= back {
            if left == back {
                self.data[left] += (1.0 - frac) * arg_value;
            }
            return;
        }
        self.data[left + 1] += frac * arg_value;
        self.data[left] += (1.0 - frac) * arg_value;
    }

    /// The interpolated derivative at outside position `arg_pos`.
    ///
    /// The derivative of the piecewise linear curve, itself interpolated
    /// linearly between the midpoints of neighbouring sample intervals, so it
    /// is continuous. The source asks its users to report that they use this.
    ///
    /// Native detail: with fewer than two samples, or at a NaN, infinite or
    /// out-of-`isize` position, a sample the source would read out of range
    /// counts as zero and a non-finite position yields NaN (for NaN) or zero.
    /// The source reads `data[1]` of a one-sample container and `data[-1]` when
    /// `left` rounds to the size of an empty one.
    pub fn derivative(&self, arg_pos: f64) -> f64 {
        let pos = self.key2index(arg_pos);
        if pos.is_nan() {
            return f64::NAN;
        }
        let sample = |i: isize| -> f64 {
            usize::try_from(i)
                .ok()
                .and_then(|i| self.data.get(i))
                .copied()
                .unwrap_or(0.0)
        };
        let size = isize::try_from(self.data.len()).unwrap_or(isize::MAX);
        // Source `int(pos + 0.5)`: truncation toward zero.
        let Some((left, _)) = split(pos + 0.5) else {
            return 0.0;
        };
        if left < 0 {
            return 0.0;
        }
        if left == 0 {
            if pos >= -0.5 {
                return (sample(1) - sample(0)) * (pos + 0.5) + sample(0) * (0.5 - pos);
            }
            return sample(0) * (pos + 1.5);
        }
        let factor = left as f64 - pos + 0.5;
        if left > size {
            return 0.0;
        }
        if left < size - 1 {
            return (sample(left) - sample(left - 1)) * factor
                + (sample(left + 1) - sample(left)) * (1.0 - factor);
        }
        if left == size - 1 {
            return (sample(left) - sample(left - 1)) * factor + (-sample(left)) * (1.0 - factor);
        }
        (-sample(left - 1)) * factor
    }

    /// The samples values are interpolated from.
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    /// Mutable access to the samples (source non-const `getData`).
    pub fn data_mut(&mut self) -> &mut Vec<f64> {
        &mut self.data
    }

    /// Replace the samples (source `setData`).
    pub fn set_data(&mut self, data: Vec<f64>) {
        self.data = data;
    }

    /// Whether there are no samples.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The transformation from outside to inside coordinates:
    /// `(pos - offset) / scale`, or `0` when the scale is zero.
    ///
    /// A NaN scale is not zero, as in the source's `if (scale_)`.
    pub fn key2index(&self, pos: f64) -> f64 {
        if self.scale != 0.0 {
            (pos - self.offset) / self.scale
        } else {
            0.0
        }
    }

    /// The transformation from inside to outside coordinates:
    /// `pos * scale + offset`.
    pub fn index2key(&self, pos: f64) -> f64 {
        pos * self.scale + self.offset
    }

    /// The scale: the distance, in outside units, between consecutive samples.
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Set the scale.
    ///
    /// As in the source, this invalidates the inside and outside reference
    /// points, which keep their old values.
    pub fn set_scale(&mut self, scale: f64) {
        self.scale = scale;
    }

    /// The offset: the outside position that corresponds to sample 0.
    pub fn offset(&self) -> f64 {
        self.offset
    }

    /// Set the offset.
    ///
    /// As in the source, this invalidates the inside and outside reference
    /// points, which keep their old values.
    pub fn set_offset(&mut self, offset: f64) {
        self.offset = offset;
    }

    /// Specify the mapping from outside to inside coordinates by `scale` (the
    /// outside distance between consecutive samples) and one pair of positions,
    /// `inside` and `outside`, that are mapped onto each other.
    ///
    /// For example, a probability density centred on zero cannot have negative
    /// indices; precompute it shifted, and map its mean (`inside`) onto the
    /// centroid of a real peak (`outside`) to model that peak by a scaled and
    /// shifted copy of the density. The offset becomes
    /// `outside - scale * inside`.
    pub fn set_mapping(&mut self, scale: f64, inside: f64, outside: f64) {
        self.scale = scale;
        self.inside = inside;
        self.outside = outside;
        self.offset = outside - scale * inside;
    }

    /// Specify the mapping by two pairs of positions that are mapped onto each
    /// other: `inside_low` onto `outside_low` and `inside_high` onto
    /// `outside_high`.
    ///
    /// A convenience form of [`Self::set_mapping`]: the scale is
    /// `(outside_high - outside_low) / (inside_high - inside_low)`, or zero
    /// when the two inside positions are equal (source four-argument
    /// `setMapping`).
    pub fn set_mapping_by_points(
        &mut self,
        inside_low: f64,
        outside_low: f64,
        inside_high: f64,
        outside_high: f64,
    ) {
        if inside_high != inside_low {
            self.set_mapping(
                (outside_high - outside_low) / (inside_high - inside_low),
                inside_low,
                outside_low,
            );
        } else {
            self.set_mapping(0.0, inside_low, outside_low);
        }
    }

    /// The inside reference point of the last [`Self::set_mapping`].
    pub fn inside_reference_point(&self) -> f64 {
        self.inside
    }

    /// The outside reference point of the last [`Self::set_mapping`].
    pub fn outside_reference_point(&self) -> f64 {
        self.outside
    }

    /// Lower boundary of the support, in outside coordinates: the position of
    /// inside index -1, or of index 0 when there is no data.
    pub fn support_min(&self) -> f64 {
        self.index2key(if self.is_empty() { 0.0 } else { -1.0 })
    }

    /// Upper boundary of the support, in outside coordinates: the position of
    /// inside index `data.len()`.
    pub fn support_max(&self) -> f64 {
        self.index2key(self.data.len() as f64)
    }
}
