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

/// Raw RT/m/z point. Public array storage is the safe mutable DPosition<2>
/// equivalent; coordinates are f64 and intensity is f32.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Peak2D {
    pub position: [f64; 2],
    pub intensity: f32,
}
/// Raw mobility/m/z point. The value does not carry a mobility unit or scan RT.
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
            pub const $dimension: usize = 0;
            pub const MZ: usize = 1;
            pub const DIMENSION: usize = 2;
            pub const fn new(first: f64, mz: f64, intensity: f32) -> Self {
                Self {
                    position: [first, mz],
                    intensity,
                }
            }
            pub const fn from_position(position: [f64; 2], intensity: f32) -> Self {
                Self {
                    position,
                    intensity,
                }
            }
            pub const fn $get(&self) -> f64 {
                self.position[0]
            }
            pub fn $set(&mut self, value: f64) {
                self.position[0] = value;
            }
            pub const fn mz(&self) -> f64 {
                self.position[1]
            }
            pub fn set_mz(&mut self, value: f64) {
                self.position[1] = value;
            }
            pub fn short_dimension_name(dimension: usize) -> Result<&'static str> {
                dimension_value([$label, "MZ"], dimension)
            }
            pub fn full_dimension_name(dimension: usize) -> Result<&'static str> {
                dimension_value([$name, "mass-to-charge"], dimension)
            }
            pub fn short_dimension_unit(dimension: usize) -> Result<&'static str> {
                dimension_value([$short_unit, "Th"], dimension)
            }
            pub fn full_dimension_unit(dimension: usize) -> Result<&'static str> {
                dimension_value([$full_unit, "Thomson"], dimension)
            }
            pub const fn $short_name() -> &'static str {
                $label
            }
            pub const fn short_dimension_name_mz() -> &'static str {
                "MZ"
            }
            pub const fn $full_name() -> &'static str {
                $name
            }
            pub const fn full_dimension_name_mz() -> &'static str {
                "mass-to-charge"
            }
            pub const fn $short_unit_fn() -> &'static str {
                $short_unit
            }
            pub const fn short_dimension_unit_mz() -> &'static str {
                "Th"
            }
            pub const fn $full_unit_fn() -> &'static str {
                $full_unit
            }
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

/// Annotated RT/m/z point. Own value equality includes metadata and unique ID;
/// metadata follows the existing native exact-value/unit contract.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct RichPeak2D {
    pub peak: Peak2D,
    pub metadata: MetaInfo,
    /// Zero is the source invalid/unassigned value. Generation is separate.
    pub unique_id: u64,
}
impl RichPeak2D {
    pub fn new(rt: f64, mz: f64, intensity: f32) -> Self {
        Peak2D::new(rt, mz, intensity).into()
    }
    pub fn from_position(position: [f64; 2], intensity: f32) -> Self {
        Peak2D::from_position(position, intensity).into()
    }
    /// Assign from a plain point, clearing current metadata and ID. Return old
    /// ownership so this constant-time operation does not destroy its payload.
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
