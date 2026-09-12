// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Member-review closures for `KERNEL/FeatureHandle.h` and
//! `KERNEL/BinnedSpectrum.h`.
//!
//! The member-by-member review recorded in `docs/FEATURE_SUPPORT.md` and
//! `docs/COMPARISON_SUPPORT.md` found four public source members without a
//! Rust counterpart. They are implemented here rather than in the types' home
//! files (`kernel/features.rs`, `comparison.rs`) because those files were
//! frozen for the review work package; Rust permits inherent and trait
//! implementations for a crate-local type from any module of the same crate,
//! so callers see no difference. Nothing here changes existing behaviour.

use super::features::FeatureHandle;
use super::peak2d::Peak2D;
use crate::comparison::BinnedSpectrum;
use crate::concept::HasUniqueId;
use std::fmt;
use std::hash::{Hash, Hasher};

impl FeatureHandle {
    /// Handle for a plain point in map `map_index`, identified by
    /// `element_index`.
    ///
    /// Source `FeatureHandle(UInt64 map_index, const Peak2D& point, UInt64
    /// element_index)`: the point supplies RT, m/z and intensity, the element
    /// index becomes the unique ID, and charge and width start at zero. The
    /// value is not validated; call [`FeatureHandle::validate`] before storing
    /// it in a consensus feature, which validates on insertion anyway.
    pub fn from_peak(map_index: u64, point: Peak2D, element_index: u64) -> Self {
        Self {
            map_index,
            unique_id: element_index,
            rt: point.rt(),
            mz: point.mz(),
            intensity: point.intensity,
            charge: 0,
            width: 0.0,
        }
    }
}

/// The inherited `UniqueIdInterface` surface over the public `unique_id`
/// field: validity checks, clear, swap, suffix-text parsing and caller-owned
/// generator assignment.
impl HasUniqueId for FeatureHandle {
    fn unique_id(&self) -> u64 {
        self.unique_id
    }
    fn unique_id_mut(&mut self) -> &mut u64 {
        &mut self.unique_id
    }
}

/// Source `operator<<(std::ostream&, const FeatureHandle&)` layout: a banner
/// line followed by RT, m/z, intensity, map index and element ID, one per
/// line, each terminated by a newline. Charge and width are not printed, as
/// in the source. Numbers use Rust formatting rather than the C++ stream's
/// six-significant-digit default, matching every other kernel `Display`.
impl fmt::Display for FeatureHandle {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "---------- FeatureHandle -----------------")?;
        writeln!(f, "RT: {}", self.rt)?;
        writeln!(f, "m/z: {}", self.mz)?;
        writeln!(f, "Intensity: {}", self.intensity)?;
        writeln!(f, "Map Index: {}", self.map_index)?;
        writeln!(f, "Element Id: {}", self.unique_id)
    }
}

/// Source `std::hash<FeatureHandle>` hashes RT, m/z, intensity, unique ID,
/// map index, charge and width. This feeds the same seven values to the
/// caller's `Hasher`, normalising both zero signs of each float as the source
/// `hash_float` does, so `PartialEq`-equal handles hash equally. No digest
/// value is shared with the C++ FNV-1a combination.
impl Hash for FeatureHandle {
    fn hash<H: Hasher>(&self, state: &mut H) {
        float_bits(self.rt).hash(state);
        float_bits(self.mz).hash(state);
        float_bits(f64::from(self.intensity)).hash(state);
        self.unique_id.hash(state);
        self.map_index.hash(state);
        self.charge.hash(state);
        float_bits(f64::from(self.width)).hash(state);
    }
}

fn float_bits(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}

/// Recommended bin layouts from doi:10.1007/s13361-015-1179-x, as declared on
/// the source class.
///
/// The source `@todo` notes that weighted intensity spread (half intensity in
/// each flanking bin for high-resolution sum scores) is not implemented; this
/// port does not implement it either, so `spread` remains an integer bin count.
impl BinnedSpectrum {
    /// Source `DEFAULT_BIN_WIDTH_LOWRES`; also the [`crate::comparison::BinConfig`]
    /// default `size`.
    pub const DEFAULT_BIN_WIDTH_LOWRES: f32 = 1.0005;
    /// Source `DEFAULT_BIN_WIDTH_HIRES`.
    pub const DEFAULT_BIN_WIDTH_HIRES: f32 = 0.02;
    /// Source `DEFAULT_BIN_OFFSET_HIRES`.
    pub const DEFAULT_BIN_OFFSET_HIRES: f32 = 0.0;
    /// Source `DEFAULT_BIN_OFFSET_LOWRES`; also the [`crate::comparison::BinConfig`]
    /// default `offset`.
    pub const DEFAULT_BIN_OFFSET_LOWRES: f32 = 0.4;
}
