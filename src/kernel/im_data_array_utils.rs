// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Interpreting ion-mobility float data array names, independently of frame
//! conversion.
//!
//! Ports `src/openms/include/OpenMS/IONMOBILITY/IMDataArrayUtils.h` and
//! `src/openms/source/IONMOBILITY/IMDataArrayUtils.cpp` (core SDK `bc9cc12`).
//! See `docs/SPECTRUM_MOBILITY_SUPPORT.md`.
//!
//! An ion-mobility frame annotates every peak through a float data array whose
//! *name* says it holds ion mobility, and in which unit: a child term of
//! `MS:1002893 ! ion mobility array`, or one of the vendor `UserParam` names.
//! [`ImDataArrayUtils::set_im_unit`] writes such a name for a unit and
//! [`ImDataArrayUtils::im_unit`] reads the unit back. `IMDataConverter` keeps
//! both as compatibility entry points, which
//! `crate::kernel::im_data_converter::ImDataConverter` forwards here too.
//!
//! The source resolves names against the PSI-MS vocabulary at run time; this
//! port pins the nine child terms and their units from `resources/cv/psi-ms.obo`
//! (see `crate::kernel::spectrum_mobility`), so no controlled-vocabulary file is
//! needed at run time.

use super::DataArray;
use super::spectrum_mobility::array_unit;
use crate::metadata::DriftTimeUnit;
use crate::{Error, Result};

/// `MS:1002816 ! mean ion mobility array`, the name
/// [`ImDataArrayUtils::set_im_unit`] writes for milliseconds.
pub const MEAN_ION_MOBILITY_ARRAY: &str = "mean ion mobility array";

/// `MS:1003008 ! raw inverse reduced ion mobility array`, the name
/// [`ImDataArrayUtils::set_im_unit`] writes for `1/K0` (volt-second per square
/// centimetre).
pub const RAW_INVERSE_REDUCED_ION_MOBILITY_ARRAY: &str = "raw inverse reduced ion mobility array";

/// The source class `IMDataArrayUtils`, whose two members are static.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ImDataArrayUtils;

impl ImDataArrayUtils {
    /// Name `array` after the ion-mobility unit `unit` (`setIMUnit`).
    ///
    /// | unit | name written |
    /// |---|---|
    /// | [`DriftTimeUnit::Millisecond`] | `mean ion mobility array` (`MS:1002816`) |
    /// | [`DriftTimeUnit::InverseReducedMobility`] | `raw inverse reduced ion mobility array` (`MS:1003008`) |
    ///
    /// The mapping is not exact — the name cannot say whether the values are
    /// raw or binned — but it lets [`ImDataArrayUtils::im_unit`] recover the
    /// unit. Only the name changes; the values and the rest of the array are
    /// untouched.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other unit ([`DriftTimeUnit::None`],
    /// FAIMS compensation voltage, collision cross section), where the source
    /// throws `Exception::InvalidValue` ("Unit is not a valid IM unit for float
    /// data arrays"): the PSI-MS vocabulary has no ion-mobility array term for
    /// them. The array is left unchanged.
    pub fn set_im_unit(array: &mut DataArray<f32>, unit: DriftTimeUnit) -> Result<()> {
        let name = match unit {
            DriftTimeUnit::Millisecond => MEAN_ION_MOBILITY_ARRAY,
            DriftTimeUnit::InverseReducedMobility => RAW_INVERSE_REDUCED_ION_MOBILITY_ARRAY,
            other => {
                return Err(Error::InvalidValue(format!(
                    "Unit is not a valid IM unit for float data arrays: '{}'",
                    other.name()
                )));
            }
        };
        array.name = name.to_string();
        Ok(())
    }

    /// The ion-mobility unit `array`'s name denotes, or `None` when the array
    /// does not hold ion mobility (`getIMUnit`).
    ///
    /// The name is looked up as the source does, first match wins:
    ///
    /// 1. an exact PSI-MS child of `MS:1002893 ! ion mobility array`, whose
    ///    declared unit decides — `MS:1002814` (volt-second per square
    ///    centimetre) gives [`DriftTimeUnit::InverseReducedMobility`],
    ///    `UO:0000028` [`DriftTimeUnit::Millisecond`];
    /// 2. a name starting with one of the two inverse-reduced `UserParam`
    ///    names gives [`DriftTimeUnit::InverseReducedMobility`];
    /// 3. a name starting with the generic `Ion Mobility` `UserParam` name
    ///    gives `1/K0` when it embeds `MS:1002815` or `MS:1003006`, collision
    ///    cross section when it embeds `MS:1002954`, and milliseconds
    ///    otherwise.
    ///
    /// Prefix matching is case-sensitive. The source returns `bool` and writes
    /// the unit through an out-parameter left unchanged for a non-IM array;
    /// `Option` carries both.
    ///
    /// # Notes
    ///
    /// The source also answers `true` with [`DriftTimeUnit::None`] for a CV
    /// child term that declares none of the three units. At the pinned CV every
    /// one of the nine children declares one, so that branch is unreachable and
    /// has no counterpart; `docs/SPECTRUM_MOBILITY_SUPPORT.md` records this.
    pub fn im_unit(array: &DataArray<f32>) -> Option<DriftTimeUnit> {
        array_unit(&array.name)
    }
}
