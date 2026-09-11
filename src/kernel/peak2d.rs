// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use crate::concept::HasUniqueId;
use crate::metadata::MetaInfo;
use crate::{Error, Result};
use std::{
    fmt,
    hash::{Hash, Hasher},
    ops::{Deref, DerefMut},
};

/// A two-dimensional raw data point or peak, over retention time and m/z.
///
/// Intended for continuous data as well as peak data. To annotate a single peak
/// with metadata, use [`RichPeak2D`] instead.
///
/// The public `position` array is the safe mutable equivalent of the source
/// `DPosition<2>`: coordinates are `f64` and intensity is `f32`, as in source.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Peak2D {
    pub position: [f64; 2],
    pub intensity: f32,
}
/// A two-dimensional raw data point or peak, over ion mobility and m/z.
///
/// The value carries neither a mobility unit nor the scan retention time; both
/// belong to the enclosing record, as in source.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MobilityPeak2D {
    pub position: [f64; 2],
    pub intensity: f32,
}

macro_rules! point_surface {
    ($ty:ident, $dimension:ident, $get:ident, $set:ident, $label:literal,
     $name:literal, $short_unit:literal, $full_unit:literal,
     $short_name:ident, $full_name:ident, $short_unit_fn:ident, $full_unit_fn:ident) => {
        impl $ty {
            /// Dimension index of the first coordinate.
            pub const $dimension: usize = 0;
            /// Dimension index of the mass-to-charge coordinate.
            pub const MZ: usize = 1;
            /// The number of dimensions.
            pub const DIMENSION: usize = 2;
            /// A point at the given coordinates and intensity.
            pub const fn new(first: f64, mz: f64, intensity: f32) -> Self {
                Self {
                    position: [first, mz],
                    intensity,
                }
            }
            /// A point from a position array ordered as the dimension indices.
            pub const fn from_position(position: [f64; 2], intensity: f32) -> Self {
                Self {
                    position,
                    intensity,
                }
            }
            /// The first coordinate.
            pub const fn $get(&self) -> f64 {
                self.position[0]
            }
            /// Set the first coordinate.
            pub fn $set(&mut self, value: f64) {
                self.position[0] = value;
            }
            /// The mass-to-charge coordinate.
            pub const fn mz(&self) -> f64 {
                self.position[1]
            }
            /// Set the mass-to-charge coordinate.
            pub fn set_mz(&mut self, value: f64) {
                self.position[1] = value;
            }
            /// The abbreviated name of a dimension.
            ///
            /// # Errors
            ///
            /// Returns [`Error::InvalidValue`] when `dimension` is not a valid
            /// index. Source indexes a fixed array without checking.
            pub fn short_dimension_name(dimension: usize) -> Result<&'static str> {
                dimension_value([$label, "MZ"], dimension)
            }
            /// The self-explanatory name of a dimension.
            ///
            /// # Errors
            ///
            /// Returns [`Error::InvalidValue`] when `dimension` is not a valid
            /// index. Source indexes a fixed array without checking.
            pub fn full_dimension_name(dimension: usize) -> Result<&'static str> {
                dimension_value([$name, "mass-to-charge"], dimension)
            }
            /// The abbreviated unit of measurement of a dimension.
            ///
            /// # Errors
            ///
            /// Returns [`Error::InvalidValue`] when `dimension` is not a valid
            /// index. Source indexes a fixed array without checking.
            pub fn short_dimension_unit(dimension: usize) -> Result<&'static str> {
                dimension_value([$short_unit, "Th"], dimension)
            }
            /// The self-explanatory unit of measurement of a dimension.
            ///
            /// # Errors
            ///
            /// Returns [`Error::InvalidValue`] when `dimension` is not a valid
            /// index. Source indexes a fixed array without checking.
            pub fn full_dimension_unit(dimension: usize) -> Result<&'static str> {
                dimension_value([$full_unit, "Thomson"], dimension)
            }
            /// The abbreviated name of the first dimension.
            pub const fn $short_name() -> &'static str {
                $label
            }
            /// The abbreviated name of the m/z dimension.
            pub const fn short_dimension_name_mz() -> &'static str {
                "MZ"
            }
            /// The self-explanatory name of the first dimension.
            pub const fn $full_name() -> &'static str {
                $name
            }
            /// The self-explanatory name of the m/z dimension.
            pub const fn full_dimension_name_mz() -> &'static str {
                "mass-to-charge"
            }
            /// The abbreviated unit of the first dimension.
            pub const fn $short_unit_fn() -> &'static str {
                $short_unit
            }
            /// The abbreviated unit of the m/z dimension.
            pub const fn short_dimension_unit_mz() -> &'static str {
                "Th"
            }
            /// The self-explanatory unit of the first dimension.
            pub const fn $full_unit_fn() -> &'static str {
                $full_unit
            }
            /// The self-explanatory unit of the m/z dimension.
            pub const fn full_dimension_unit_mz() -> &'static str {
                "Thomson"
            }
        }
        impl Hash for $ty {
            fn hash<H: Hasher>(&self, state: &mut H) {
                for value in self.position {
                    (if value == 0. { 0 } else { value.to_bits() }).hash(state);
                }
                (if self.intensity == 0. {
                    0
                } else {
                    self.intensity.to_bits()
                })
                .hash(state);
            }
        }
        impl fmt::Display for $ty {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                match f.precision() {
                    Some(p) => write!(
                        f,
                        concat!($label, ": {:.*} MZ: {:.*} INT: {:.*}"),
                        p, self.position[0], p, self.position[1], p, self.intensity
                    ),
                    None => write!(
                        f,
                        concat!($label, ": {} MZ: {} INT: {}"),
                        self.position[0], self.position[1], self.intensity
                    ),
                }
            }
        }
    };
}
point_surface!(
    Peak2D,
    RT,
    rt,
    set_rt,
    "RT",
    "retention time",
    "sec",
    "Seconds",
    short_dimension_name_rt,
    full_dimension_name_rt,
    short_dimension_unit_rt,
    full_dimension_unit_rt
);
point_surface!(
    MobilityPeak2D,
    IM,
    mobility,
    set_mobility,
    "IM",
    "ion mobility",
    "?",
    "?",
    short_dimension_name_im,
    full_dimension_name_im,
    short_dimension_unit_im,
    full_dimension_unit_im
);

fn dimension_value(values: [&'static str; 2], dimension: usize) -> Result<&'static str> {
    values
        .get(dimension)
        .copied()
        .ok_or_else(|| Error::InvalidValue("2D peak dimension must be 0 or 1".into()))
}

/// A two-dimensional raw data point or peak with meta information.
///
/// Intended for continuous data as well as peak data. When single peaks need no
/// metadata, use [`Peak2D`] instead.
///
/// Value equality includes the metadata and the unique ID, and the metadata
/// follows the crate's exact-value and unit contract.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct RichPeak2D {
    pub peak: Peak2D,
    pub metadata: MetaInfo,
    /// Zero is the source invalid/unassigned value. Generation is separate.
    pub unique_id: u64,
}
impl RichPeak2D {
    /// An annotated point at the given retention time, m/z and intensity.
    pub fn new(rt: f64, mz: f64, intensity: f32) -> Self {
        Peak2D::new(rt, mz, intensity).into()
    }
    /// An annotated point from a `[rt, mz]` position array.
    pub fn from_position(position: [f64; 2], intensity: f32) -> Self {
        Peak2D::from_position(position, intensity).into()
    }
    /// Replace the coordinates from a plain point, clearing the metadata and
    /// unique ID, and return the previous value.
    ///
    /// Source assignment from a `Peak2D` discards the annotations in place;
    /// returning the old value hands that payload back to the caller instead of
    /// destroying it, which keeps the operation constant time.
    pub fn replace_from_peak(&mut self, peak: Peak2D) -> Self {
        std::mem::replace(self, peak.into())
    }
}
impl HasUniqueId for RichPeak2D {
    fn unique_id(&self) -> u64 {
        self.unique_id
    }
    fn unique_id_mut(&mut self) -> &mut u64 {
        &mut self.unique_id
    }
}
impl From<Peak2D> for RichPeak2D {
    fn from(peak: Peak2D) -> Self {
        Self {
            peak,
            metadata: MetaInfo::new(),
            unique_id: 0,
        }
    }
}
impl From<&Peak2D> for RichPeak2D {
    fn from(peak: &Peak2D) -> Self {
        (*peak).into()
    }
}
impl Deref for RichPeak2D {
    type Target = Peak2D;
    fn deref(&self) -> &Peak2D {
        &self.peak
    }
}
impl DerefMut for RichPeak2D {
    fn deref_mut(&mut self) -> &mut Peak2D {
        &mut self.peak
    }
}
impl fmt::Display for RichPeak2D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.peak.fmt(f)
    }
}
