// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ion-mobility type vocabulary and the conversions over it.
//!
//! Ports `IONMOBILITY/IMTypes.h` together with `IMTypes.cpp` and
//! `IMTypesExperiment.cpp`: the three enum name tables, the six string
//! conversions, the unset drift-time sentinel, the Mason-Schamp collision
//! cross-section conversions and both `determineIMFormat` overloads. See
//! `docs/IM_TYPES_SUPPORT.md`.
//!
//! The three enums themselves were ported earlier and live beside this module:
//! [`DriftTimeUnit`](crate::metadata::DriftTimeUnit) is the source
//! `DriftTimeUnit`, [`IonMobilityFormat`](crate::metadata::IonMobilityFormat)
//! is `IMFormat` and
//! [`IonMobilityPeakType`](crate::metadata::IonMobilityPeakType) is
//! `IMPeakType`. None of them carries the source's `SIZE_OF_*` sentinel value,
//! which C++ needs only to size the name arrays and to terminate its linear
//! searches; every `@throws Exception::InvalidValue if @p value is SIZE_OF_*`
//! clause of the header therefore describes a state this port cannot represent,
//! and those conversions are infallible here.
//!
//! The source writes two log records this port does not emit, because no kernel
//! or metadata module in the crate is wired to
//! [`LogStream`](crate::concept::log_stream::LogStream): a debug record when a
//! spectrum carries both a single drift time and an ion-mobility array, and a
//! warning when a spectrum has a drift time but no unit. Both conditions stay
//! observable at the call site — `spectrum.contains_im_data() &&
//! spectrum.has_drift_time()` and `spectrum.drift_time_unit ==
//! DriftTimeUnit::None` — so nothing about the determination is hidden.

use super::{DriftTimeUnit, IonMobilityFormat, IonMobilityPeakType};
use crate::kernel::{MSExperiment, MSSpectrum};
use crate::{Error, Result};
use std::collections::BTreeSet;

/// Names of the drift time units, the source `NamesOfDriftTimeUnit`
/// (`IMTypes.cpp:23`), in enum order and usable as an axis annotation.
///
/// The source array is sized by `SIZE_OF_DRIFTTIMEUNIT` and therefore holds one
/// entry per real unit; the sentinel itself has no name. Each entry equals
/// [`DriftTimeUnit::name`](crate::metadata::DriftTimeUnit::name) of the
/// corresponding variant, which `tests/im_types.rs` asserts so the two tables
/// cannot drift apart.
pub const NAMES_OF_DRIFT_TIME_UNIT: [&str; 5] = ["<NONE>", "ms", "1/K0", "FAIMS_CV", "CCS"];

/// Names of the ion-mobility formats, the source `NamesOfIMFormat`
/// (`IMTypes.cpp:24`), in enum order.
pub const NAMES_OF_IM_FORMAT: [&str; 4] = ["none", "im_peak", "im_spectrum", "unknown"];

/// Names of the ion-mobility peak types, the source `NamesOfIMPeakType`
/// (`IMTypes.cpp:69`), in enum order.
pub const NAMES_OF_IM_PEAK_TYPE: [&str; 3] = ["im_profile", "im_centroided", "unknown"];

/// Convert an entry of [`NAMES_OF_DRIFT_TIME_UNIT`] to its unit.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `name` is not one of the names, as the
/// source `toDriftTimeUnit` throws `Exception::InvalidValue` with the offending
/// string as the value. Matching is exact and case-sensitive, as the source
/// `std::find` over the array.
pub fn to_drift_time_unit(name: &str) -> Result<DriftTimeUnit> {
    DriftTimeUnit::ALL
        .iter()
        .copied()
        .find(|unit| unit.name() == name)
        .ok_or_else(|| Error::InvalidValue(format!("unknown drift time unit '{name}'")))
}

/// The name of a drift time unit, the source `driftTimeUnitToString`.
///
/// Infallible: the source's `@throws` clause fires only for
/// `SIZE_OF_DRIFTTIMEUNIT`, which [`DriftTimeUnit`](crate::metadata::DriftTimeUnit)
/// does not have.
pub fn drift_time_unit_to_string(unit: DriftTimeUnit) -> &'static str {
    unit.name()
}

/// Convert an entry of [`NAMES_OF_IM_FORMAT`] to its format.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `name` is not one of the names, as the
/// source `toIMFormat`.
pub fn to_im_format(name: &str) -> Result<IonMobilityFormat> {
    IonMobilityFormat::ALL
        .iter()
        .copied()
        .find(|format| format.name() == name)
        .ok_or_else(|| Error::InvalidValue(format!("unknown ion mobility format '{name}'")))
}

/// The name of an ion-mobility format, the source `imFormatToString`.
///
/// Infallible, for the reason given on [`drift_time_unit_to_string`].
pub fn im_format_to_string(format: IonMobilityFormat) -> &'static str {
    format.name()
}

/// Convert an entry of [`NAMES_OF_IM_PEAK_TYPE`] to its peak type.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `name` is not one of the names, as the
/// source `toIMPeakType`.
pub fn to_im_peak_type(name: &str) -> Result<IonMobilityPeakType> {
    IonMobilityPeakType::ALL
        .iter()
        .copied()
        .find(|peak_type| peak_type.name() == name)
        .ok_or_else(|| Error::InvalidValue(format!("unknown ion mobility peak type '{name}'")))
}

/// The name of an ion-mobility peak type, the source `imPeakTypeToString`.
///
/// Infallible: the source guards `>= SIZE_OF_IMPEAKTYPE`, a value
/// [`IonMobilityPeakType`](crate::metadata::IonMobilityPeakType) cannot hold.
pub fn im_peak_type_to_string(peak_type: IonMobilityPeakType) -> &'static str {
    peak_type.name()
}

/// The Bruker Mason-Schamp calibration constant relating `1/K0` in V*s/cm^2 to
/// a collision cross section in Angstrom^2 for an N2 drift gas
/// (`IMTypes.cpp:112`). Private in the source's anonymous namespace, so it is
/// private here too.
const MASON_SCHAMP_CONSTANT: f64 = 1059.62245;

/// Ion-gas reduced mass in Da, with the ion mass approximated as `mz * |charge|`
/// (source `reducedMass_`, `IMTypes.cpp:115-119`).
fn reduced_mass(mz: f64, charge: i32, buffer_gas_mass: f64) -> f64 {
    let ion_mass = mz * f64::from(charge.unsigned_abs());
    (ion_mass * buffer_gas_mass) / (ion_mass + buffer_gas_mass)
}

/// Reject the shared argument domain of the two cross-section conversions.
///
/// The source tests `value <= 0.0 || mz <= 0.0 || charge == 0`; a NaN fails
/// every one of those comparisons and is accepted, which is why this port also
/// requires finiteness.
fn check_ccs_arguments(
    what: &str,
    value: f64,
    mz: f64,
    charge: i32,
    buffer_gas_mass: f64,
) -> Result<()> {
    let positive = |name: &str, number: f64| -> Result<()> {
        if number.is_finite() && number > 0.0 {
            Ok(())
        } else {
            Err(Error::InvalidValue(format!(
                "{what} requires a finite positive {name}, got {number}"
            )))
        }
    };
    positive(what, value)?;
    positive("m/z", mz)?;
    positive("buffer gas mass", buffer_gas_mass)?;
    if charge == 0 {
        return Err(Error::InvalidValue(format!(
            "{what} requires a nonzero charge"
        )));
    }
    Ok(())
}

/// Static ion-mobility helpers, the source `class IMTypes`.
///
/// The source class has no state and is never instantiated for its own sake;
/// its class test constructs and deletes one, so this is a unit struct with a
/// [`Default`] implementation and the compiler-generated drop.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ImTypes;

impl ImTypes {
    /// The drift time of a spectrum that is not an ion-mobility spectrum, the
    /// source `IMTypes::DRIFTTIME_NOT_SET` (`IMTypes.h:93`).
    ///
    /// [`MSSpectrum::drift_time`](crate::kernel::MSSpectrum::drift_time) stores
    /// this sentinel for an unset drift time, exactly as the source does;
    /// `MSSpectrum::drift_time_if_set` is the checked accessor over it.
    pub const DRIFTTIME_NOT_SET: f64 = -1.0;

    /// Mass of N2, the default buffer gas of the cross-section conversions, in
    /// Da (source `IMTypes::N2_BUFFER_GAS_MASS`).
    ///
    /// The source uses the rounded 28.0 rather than 28.006148 deliberately, to
    /// match the calibration constant that was validated against the
    /// Bruker/alphatims and MaxQuant cross-section values.
    pub const N2_BUFFER_GAS_MASS: f64 = 28.0;

    /// Spectra one [`determine_im_format_for_ms_level`] call may scan.
    ///
    /// Native ceiling; the source scans an experiment of any size. The value
    /// matches `RangeManager::MAX_ITEMS` and `MSSpectrum::MAX_MOBILITY_ITEMS`,
    /// so an experiment whose ranges can be computed can also be asked for its
    /// ion-mobility format.
    ///
    /// [`determine_im_format_for_ms_level`]: ImTypes::determine_im_format_for_ms_level
    pub const MAX_SPECTRA: usize = 100_000_000;

    /// The ion-mobility format a spectrum's data shows, the source
    /// `IMTypes::determineIMFormat(const MSSpectrum&)`
    /// (`IMTypesExperiment.cpp:53-81`).
    ///
    /// An ion-mobility float data array yields
    /// [`IonMobilityFormat::PerPeak`](crate::metadata::IonMobilityFormat::PerPeak),
    /// otherwise a set drift time yields
    /// [`IonMobilityFormat::PerSpectrum`](crate::metadata::IonMobilityFormat::PerSpectrum),
    /// otherwise [`IonMobilityFormat::None`](crate::metadata::IonMobilityFormat::None).
    /// A spectrum carrying both is `PerPeak`: the source calls that combination
    /// valid and experimental, and logs a debug record rather than failing.
    ///
    /// This is only the data half of the source function. The source first
    /// returns the spectrum's *stored* format whenever it is not
    /// [`IonMobilityFormat::Unknown`](crate::metadata::IonMobilityFormat::Unknown);
    /// the Rust [`MSSpectrum`](crate::kernel::MSSpectrum) flattens
    /// `SpectrumSettings` and carries no stored format, so a caller holding one
    /// — on a [`SpectrumSettings`](crate::metadata::SpectrumSettings), via
    /// [`SpectrumSettings::im_format`](crate::metadata::SpectrumSettings::im_format)
    /// — passes it to
    /// [`determine_im_format_with_stored`](ImTypes::determine_im_format_with_stored)
    /// instead, which is the complete source function.
    ///
    /// The header declares `@throws Exception::InvalidValue if IM values are
    /// annotated as single drift time and float array`, but the implementation
    /// throws nothing on that input and `IMTypes_test.cpp:183-186` asserts
    /// `IM_PEAK` for it. The implementation and the test win; the header comment
    /// is recorded as a source defect in this port's support document. Because
    /// no input can fail, this returns the format directly rather than a
    /// `Result`.
    ///
    /// Detecting the array reuses
    /// [`MSSpectrum::contains_im_data`](crate::kernel::MSSpectrum::contains_im_data),
    /// which is the source `containsIMData` and inspects only float-array names,
    /// so the cost is one pass over the spectrum's data-array headers and no
    /// allocation.
    pub fn determine_im_format(spectrum: &MSSpectrum) -> IonMobilityFormat {
        if spectrum.contains_im_data() {
            IonMobilityFormat::PerPeak
        } else if spectrum.has_drift_time() {
            IonMobilityFormat::PerSpectrum
        } else {
            IonMobilityFormat::None
        }
    }

    /// The complete source `determineIMFormat(const MSSpectrum&)`, including the
    /// stored-format short circuit that the Rust spectrum cannot supply itself.
    ///
    /// `stored` is the spectrum's `getIMFormat()`: any value other than
    /// [`IonMobilityFormat::Unknown`](crate::metadata::IonMobilityFormat::Unknown)
    /// is returned unchanged, without looking at the data at all — so a stored
    /// [`IonMobilityFormat::None`](crate::metadata::IonMobilityFormat::None)
    /// suppresses detection on a spectrum that does carry ion mobility. Pass
    /// `Unknown`, the source default, to determine the format from the data, as
    /// [`determine_im_format`](ImTypes::determine_im_format) does.
    pub fn determine_im_format_with_stored(
        stored: IonMobilityFormat,
        spectrum: &MSSpectrum,
    ) -> IonMobilityFormat {
        if stored != IonMobilityFormat::Unknown {
            return stored;
        }
        Self::determine_im_format(spectrum)
    }

    /// The single ion-mobility format shared by the spectra of one MS level, the
    /// source `IMTypes::determineIMFormat(const MSExperiment&, int)`
    /// (`IMTypesExperiment.cpp:21-50`).
    ///
    /// Spectra of other MS levels are skipped, and
    /// [`IonMobilityFormat::None`](crate::metadata::IonMobilityFormat::None) is
    /// erased from the collected set before it is judged, so a level whose
    /// spectra are a mix of ion-mobility spectra and plain ones takes the format
    /// of the ion-mobility ones. A level with no spectra, or none carrying ion
    /// mobility, is `None`.
    ///
    /// Each spectrum is classified by
    /// [`determine_im_format`](ImTypes::determine_im_format), the data half of
    /// the per-spectrum function; the source classifies it by the full function,
    /// which would consult each spectrum's stored format first. The port's
    /// spectra carry no stored format, so the two agree for every spectrum whose
    /// stored format is the `UNKNOWN` default.
    ///
    /// A negative `ms_level` matches nothing: the source compares with
    /// `std::cmp_not_equal` against an unsigned level, so the sign is preserved
    /// rather than wrapped.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when spectra of `ms_level` disagree,
    /// naming the number of distinct formats as the source does, and when the
    /// experiment holds more than [`ImTypes::MAX_SPECTRA`] spectra. The ceiling
    /// is checked before anything is collected, and nothing is mutated in any
    /// case.
    pub fn determine_im_format_for_ms_level(
        experiment: &MSExperiment,
        ms_level: i32,
    ) -> Result<IonMobilityFormat> {
        if experiment.spectra.len() > Self::MAX_SPECTRA {
            return Err(Error::InvalidValue(format!(
                "experiment holds {} spectra, more than the {} this scan allows",
                experiment.spectra.len(),
                Self::MAX_SPECTRA
            )));
        }
        let level = u32::try_from(ms_level).ok();
        let mut seen: BTreeSet<IonMobilityFormat> = BTreeSet::new();
        for spectrum in &experiment.spectra {
            if level != Some(spectrum.ms_level) {
                continue;
            }
            seen.insert(Self::determine_im_format(spectrum));
        }
        seen.remove(&IonMobilityFormat::None);
        let mut formats = seen.into_iter();
        let Some(format) = formats.next() else {
            return Ok(IonMobilityFormat::None);
        };
        if formats.next().is_some() {
            return Err(Error::InvalidValue(format!(
                "experiment contains MS{ms_level} spectra with different ion mobility formats; handle per spectrum"
            )));
        }
        // The source guards the single-format case against anything but IM_PEAK
        // and IM_SPECTRUM and throws "subfunction returned invalid value(s)".
        // The guard is unreachable in C++ and here: the per-spectrum function
        // returns only NONE, IM_PEAK or IM_SPECTRUM and NONE was just erased.
        // It is reproduced so that a future change to the classification cannot
        // silently return an unusable format.
        if !matches!(
            format,
            IonMobilityFormat::PerPeak | IonMobilityFormat::PerSpectrum
        ) {
            return Err(Error::InvalidValue(format!(
                "ion mobility format determination returned {format}"
            )));
        }
        Ok(format)
    }

    /// Convert a reduced inverse ion mobility to a collision cross section for
    /// an N2 drift gas, the source `IMTypes::oneOverK0ToCCS` with its default
    /// `buffer_gas_mass`.
    ///
    /// See
    /// [`one_over_k0_to_ccs_with_buffer_gas`](ImTypes::one_over_k0_to_ccs_with_buffer_gas)
    /// for the relation, the arguments and the errors; this passes
    /// [`ImTypes::N2_BUFFER_GAS_MASS`].
    pub fn one_over_k0_to_ccs(one_over_k0: f64, mz: f64, charge: i32) -> Result<f64> {
        Self::one_over_k0_to_ccs_with_buffer_gas(one_over_k0, mz, charge, Self::N2_BUFFER_GAS_MASS)
    }

    /// Convert a reduced inverse ion mobility to a collision cross section via
    /// the Mason-Schamp relation, the source `IMTypes::oneOverK0ToCCS`.
    ///
    /// The conventional single-temperature form is used,
    ///
    /// ```text
    /// CCS = (C * |charge| / sqrt(mu)) * (1/K0)
    /// ```
    ///
    /// with the Bruker calibration constant `C = 1059.62245`, validated against
    /// alphatims and MaxQuant cross sections for an N2 drift gas at the usual
    /// calibration temperature, and the ion-gas reduced mass
    /// `mu = (m_ion * m_gas) / (m_ion + m_gas)` with the ion mass approximated
    /// as `m_ion = mz * |charge|`. The multiplications are written in the
    /// source's order, so the result is bit-identical to the C++ for equal
    /// inputs.
    ///
    /// # Arguments
    ///
    /// * `one_over_k0` — reduced inverse ion mobility in V*s/cm^2; must be
    ///   positive.
    /// * `mz` — precursor m/z of the ion; must be positive.
    /// * `charge` — precursor charge; the sign is ignored and it must be
    ///   nonzero.
    /// * `buffer_gas_mass` — drift-gas mass in Da; must be positive.
    ///   [`ImTypes::N2_BUFFER_GAS_MASS`] is the source default.
    ///
    /// Returns the collision cross section in square Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `one_over_k0`, `mz` or
    /// `buffer_gas_mass` is not finite and positive, or when `charge` is zero.
    /// The source checks `one_over_k0 <= 0`, `mz <= 0` and `charge == 0` only,
    /// so it accepts a NaN — every comparison against NaN is false — and
    /// returns NaN; this port rejects it, and also rejects a non-positive
    /// `buffer_gas_mass`, which the source does not examine at all even though a
    /// negative gas mass makes the reduced mass negative and the result NaN. The
    /// source's `|charge|` is `std::abs` on an `int`, which is undefined for
    /// `INT_MIN`; this uses [`i32::unsigned_abs`], so `i32::MIN` is a valid
    /// charge here.
    pub fn one_over_k0_to_ccs_with_buffer_gas(
        one_over_k0: f64,
        mz: f64,
        charge: i32,
        buffer_gas_mass: f64,
    ) -> Result<f64> {
        check_ccs_arguments("1/K0", one_over_k0, mz, charge, buffer_gas_mass)?;
        let mu = reduced_mass(mz, charge, buffer_gas_mass);
        Ok(MASON_SCHAMP_CONSTANT * f64::from(charge.unsigned_abs()) / mu.sqrt() * one_over_k0)
    }

    /// Convert a collision cross section back to a reduced inverse ion mobility
    /// for an N2 drift gas, the source `IMTypes::ccsToOneOverK0` with its
    /// default `buffer_gas_mass`.
    ///
    /// See
    /// [`ccs_to_one_over_k0_with_buffer_gas`](ImTypes::ccs_to_one_over_k0_with_buffer_gas)
    /// for the arguments and the errors; this passes
    /// [`ImTypes::N2_BUFFER_GAS_MASS`].
    pub fn ccs_to_one_over_k0(ccs: f64, mz: f64, charge: i32) -> Result<f64> {
        Self::ccs_to_one_over_k0_with_buffer_gas(ccs, mz, charge, Self::N2_BUFFER_GAS_MASS)
    }

    /// Convert a collision cross section back to a reduced inverse ion mobility,
    /// the source `IMTypes::ccsToOneOverK0` and the inverse of
    /// [`one_over_k0_to_ccs_with_buffer_gas`](ImTypes::one_over_k0_to_ccs_with_buffer_gas).
    ///
    /// # Arguments
    ///
    /// * `ccs` — collision cross section in square Angstrom; must be positive.
    /// * `mz` — precursor m/z of the ion; must be positive.
    /// * `charge` — precursor charge; the sign is ignored and it must be
    ///   nonzero.
    /// * `buffer_gas_mass` — drift-gas mass in Da; must be positive.
    ///
    /// Returns the reduced inverse ion mobility in V*s/cm^2. The round trip
    /// through both conversions is exact to within one rounding of the shared
    /// `sqrt(mu)` factor, not bit-exact, because the two expressions divide and
    /// multiply in different orders — as in the source, whose own class test
    /// compares the round trip with a tolerance.
    ///
    /// # Errors
    ///
    /// As for
    /// [`one_over_k0_to_ccs_with_buffer_gas`](ImTypes::one_over_k0_to_ccs_with_buffer_gas),
    /// with `ccs` in place of `one_over_k0`.
    pub fn ccs_to_one_over_k0_with_buffer_gas(
        ccs: f64,
        mz: f64,
        charge: i32,
        buffer_gas_mass: f64,
    ) -> Result<f64> {
        check_ccs_arguments("CCS", ccs, mz, charge, buffer_gas_mass)?;
        let mu = reduced_mass(mz, charge, buffer_gas_mass);
        Ok(ccs * mu.sqrt() / (MASON_SCHAMP_CONSTANT * f64::from(charge.unsigned_abs())))
    }
}
