// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FAIMS compensation-voltage queries, the port of `IONMOBILITY/FAIMSHelper.h`
//! and `IONMOBILITY/FAIMSHelper.cpp`.
//!
//! The source class collects the distinct compensation voltages (CVs) of an
//! experiment and filters peptide identifications by CV. Both are static
//! functions on `FaimsHelper`. The CV set is a `BTreeSet` of
//! `CompensationVoltage`, a total-ordered key that reproduces the ordering and
//! the equality of the source `std::set<double>`.
//!
//! Module-level links are plain code spans: the module's `pub mod` line carries
//! its own outer doc comment, and rustdoc resolves combined module docs in the
//! parent module, where these names are not in scope.
//!
//! This lives in `kernel` rather than in a new `ionmobility` module: the
//! other mobility operations already live here (`spectrum_mobility`,
//! `experiment_mobility`), and everything this module needs is within the
//! module edges `kernel` already has.
//!
//! See `docs/FAIMS_HELPER_SUPPORT.md` for the API mapping, the preserved source
//! conventions, the native differences and the evidence.

use super::MSExperiment;
use crate::concept::constants::user_param::FAIMS_CV;
use crate::identification::PeptideIdentification;
use crate::metadata::{DriftTimeUnit, ImTypes};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::hash::{Hash, Hasher};

/// A FAIMS compensation voltage, in volts, ordered and compared exactly as the
/// source `std::set<double>` orders and compares its elements.
///
/// # Ordering
///
/// The source set orders by `operator<` and treats two values as the same
/// element when neither is less than the other. For every value this type can
/// hold, that is:
///
/// - ascending numeric order, `-inf` first and `+inf` last, so the usual
///   negative FAIMS voltages come out as `-65 < -55 < -45`;
/// - `-0.0` and `+0.0` are one element, and a set keeps the sign of whichever
///   zero it received first, as `std::set::insert` does. [`volts`](Self::volts)
///   returns that stored sign;
/// - no tolerance: two adjacent doubles are distinct voltages.
///
/// NaN has no place in that order: `operator<` is false both ways for NaN, so
/// the source set treats it as equal to whatever element it is compared with,
/// which breaks the strict weak ordering `std::set` requires and makes the
/// result depend on insertion order. [`new`](Self::new) therefore refuses NaN,
/// and [`Ord`] is a genuine total order.
///
/// The order is implemented as [`f64::total_cmp`] on the value with negative
/// zero folded into positive zero. Plain `f64::to_bits` is *not* this order:
/// IEEE 754 bit patterns put every negative value after every positive one and
/// order negative values by magnitude, so `(-45.0f64).to_bits() <
/// (-65.0f64).to_bits()`. For FAIMS data, which is mostly negative, bit order
/// would reverse the listing FileInfo prints and the order in which
/// `IMDataConverter::splitByFAIMSCV` creates its groups.
#[derive(Clone, Copy, Debug)]
pub struct CompensationVoltage(f64);

impl CompensationVoltage {
    /// A compensation voltage of `volts`.
    ///
    /// Infinite values are accepted: `operator<` orders them, so the source
    /// set holds them without difficulty.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `volts` is NaN, which the source
    /// `std::set<double>` cannot order (see the type documentation).
    pub fn new(volts: f64) -> Result<Self> {
        if volts.is_nan() {
            return Err(Error::InvalidValue(
                "a FAIMS compensation voltage must not be NaN".into(),
            ));
        }
        Ok(Self(volts))
    }

    /// The voltage in volts, with the sign of zero it was created with.
    pub const fn volts(self) -> f64 {
        self.0
    }

    /// The comparison key: the value with `-0.0` folded into `+0.0`.
    fn key(self) -> f64 {
        if self.0 == 0.0 { 0.0 } else { self.0 }
    }
}

impl PartialEq for CompensationVoltage {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for CompensationVoltage {}

impl PartialOrd for CompensationVoltage {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CompensationVoltage {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key().total_cmp(&other.key())
    }
}

impl Hash for CompensationVoltage {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.key().to_bits().hash(state);
    }
}

impl From<CompensationVoltage> for f64 {
    fn from(voltage: CompensationVoltage) -> Self {
        voltage.volts()
    }
}

impl TryFrom<f64> for CompensationVoltage {
    type Error = Error;

    fn try_from(volts: f64) -> Result<Self> {
        Self::new(volts)
    }
}

/// The result of [`FaimsHelper::get_compensation_voltages`]: the distinct
/// voltages and the warnings the source logs while collecting them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompensationVoltages {
    /// The distinct compensation voltages in ascending order, the source return
    /// value. Empty when the experiment holds no FAIMS spectrum.
    pub voltages: BTreeSet<CompensationVoltage>,
    /// The records the source writes with `OPENMS_LOG_WARN` instead of
    /// returning them. At most one entry,
    /// [`FaimsHelper::MISSING_VOLTAGE_WARNING`], present exactly when a FAIMS
    /// spectrum carried the unset sentinel and it was removed.
    ///
    /// This kernel module is not wired to a log stream, so the caller decides
    /// where the warning goes; nothing is printed.
    pub warnings: Vec<String>,
}

impl CompensationVoltages {
    /// The voltages in volts, ascending, as the source `std::set<double>`
    /// iterates them.
    pub fn values(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.voltages.iter().map(|voltage| voltage.volts())
    }
}

/// Helper functions for FAIMS data, the source `class FAIMSHelper`.
///
/// The source documentation: *FAIMSHelper contains convenience functions to
/// deal with FAIMS compensation voltages and related data.*
///
/// The source class holds no state, and its two members are static; its class
/// test nevertheless constructs and deletes an instance
/// (`FAIMSHelper_test.cpp:32-39`). This is therefore a zero-sized unit struct:
/// [`Default`] is the constructor, and it has no drop glue, which is the
/// (virtual, empty) destructor.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FaimsHelper;

impl FaimsHelper {
    /// Spectra one [`get_compensation_voltages`](Self::get_compensation_voltages)
    /// call may scan.
    ///
    /// Native ceiling; the source scans an experiment of any size. The value is
    /// [`ImTypes::MAX_SPECTRA`], so any experiment whose ion-mobility format
    /// can be determined can also be asked for its voltages.
    pub const MAX_SPECTRA: usize = ImTypes::MAX_SPECTRA;

    /// Peptide identifications one
    /// [`filter_peptides_by_faims_cv`](Self::filter_peptides_by_faims_cv) call
    /// may examine.
    ///
    /// Native ceiling; the source filters a list of any size. It equals
    /// [`MAX_SPECTRA`](Self::MAX_SPECTRA), so identifications at most one per
    /// scannable spectrum always fit.
    pub const MAX_PEPTIDE_IDENTIFICATIONS: usize = ImTypes::MAX_SPECTRA;

    /// The source default of the `cv_tolerance` argument of
    /// `filterPeptidesByFAIMSCV` (`FAIMSHelper.h:57`), in volts.
    ///
    /// Rust has no default arguments; pass this constant for the source
    /// default.
    pub const DEFAULT_CV_TOLERANCE: f64 = 0.01;

    /// The warning the source logs when it removes the unset sentinel from the
    /// voltage set (`FAIMSHelper.cpp:52`), verbatim.
    pub const MISSING_VOLTAGE_WARNING: &'static str =
        "Warning: FAIMS compensation voltage is missing for at least one spectrum!";

    /// Get all unique FAIMS compensation voltages (CVs) that occur in an
    /// experiment, the source `getCompensationVoltages(const PeakMap&)`.
    ///
    /// The source documentation, carried across:
    ///
    /// - All spectra are scanned, and the CVs are collected from the spectra
    ///   whose drift time unit is
    ///   [`DriftTimeUnit::FaimsCompensationVoltage`]. A FAIMS spectrum anywhere
    ///   in the experiment counts, not only the first one.
    /// - If the data does not contain any FAIMS spectra, the set is empty.
    /// - The sentinel [`ImTypes::DRIFTTIME_NOT_SET`] is ignored; a warning is
    ///   logged if it was encountered. Here the warning is returned in
    ///   [`CompensationVoltages::warnings`] instead.
    ///
    /// The voltages are in volts, as the unit states. A spectrum of another
    /// unit contributes nothing, whatever its drift time, including the
    /// sentinel and NaN.
    ///
    /// The source first returns early for an experiment without spectra and
    /// then pre-scans for any FAIMS spectrum before collecting. Collecting only
    /// from FAIMS spectra yields an empty set in both of those cases anyway, so
    /// this makes one pass; the results are identical.
    ///
    /// The sentinel is removed after collection with set equality, so exactly
    /// the value `-1.0` is removed, as the source's `erase` does; `-1.0` with a
    /// FAIMS unit is never a real voltage here.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the experiment holds more than
    /// [`MAX_SPECTRA`](Self::MAX_SPECTRA) spectra, checked before anything is
    /// collected, and when a FAIMS spectrum carries a NaN drift time. The
    /// source inserts the NaN into its `std::set<double>`, whose ordering it
    /// breaks: depending on insertion order the NaN is kept, displaces nothing,
    /// or causes later voltages to be dropped. This refuses the input instead,
    /// naming the spectrum index. The experiment is only borrowed, so nothing
    /// is mutated in any case.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::kernel::faims_helper::FaimsHelper;
    /// use openms::metadata::DriftTimeUnit;
    /// use openms::{MSExperiment, MSSpectrum};
    ///
    /// let mut experiment = MSExperiment::new();
    /// for volts in [-45.0, -65.0, -45.0] {
    ///     let mut spectrum = MSSpectrum::default();
    ///     spectrum.drift_time = volts;
    ///     spectrum.drift_time_unit = DriftTimeUnit::FaimsCompensationVoltage;
    ///     experiment.spectra.push(spectrum);
    /// }
    /// let cvs = FaimsHelper::get_compensation_voltages(&experiment)?;
    /// assert_eq!(cvs.values().collect::<Vec<_>>(), [-65.0, -45.0]);
    /// assert!(cvs.warnings.is_empty());
    /// # Ok::<(), openms::Error>(())
    /// ```
    pub fn get_compensation_voltages(experiment: &MSExperiment) -> Result<CompensationVoltages> {
        if experiment.spectra.len() > Self::MAX_SPECTRA {
            return Err(Error::InvalidValue(format!(
                "experiment holds {} spectra, more than the {} this scan allows",
                experiment.spectra.len(),
                Self::MAX_SPECTRA
            )));
        }
        let mut voltages = BTreeSet::new();
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            if spectrum.drift_time_unit != DriftTimeUnit::FaimsCompensationVoltage {
                continue;
            }
            let voltage = CompensationVoltage::new(spectrum.drift_time).map_err(|_| {
                Error::InvalidValue(format!(
                    "spectrum {index} carries a NaN FAIMS compensation voltage"
                ))
            })?;
            voltages.insert(voltage);
        }
        let mut warnings = Vec::new();
        if voltages.remove(&CompensationVoltage(ImTypes::DRIFTTIME_NOT_SET)) {
            warnings.push(Self::MISSING_VOLTAGE_WARNING.to_owned());
        }
        Ok(CompensationVoltages { voltages, warnings })
    }

    /// Filter peptide identifications by FAIMS compensation voltage, the source
    /// `filterPeptidesByFAIMSCV(peptides, target_cv, cv_tolerance)`.
    ///
    /// The source documentation, carried across: the identifications are
    /// filtered to only those matching `target_cv`, and identifications without
    /// a `FAIMS_CV` annotation are included for backward compatibility.
    ///
    /// # Arguments
    ///
    /// - `peptides`: the input identifications, borrowed; the kept ones are
    ///   cloned into the result in input order.
    /// - `target_cv`: the target FAIMS compensation voltage, in volts.
    /// - `cv_tolerance`: the tolerance for the floating-point comparison, in
    ///   volts; the source default is
    ///   [`DEFAULT_CV_TOLERANCE`](Self::DEFAULT_CV_TOLERANCE) (0.01).
    ///
    /// An identification is kept when its metadata has no
    /// [`FAIMS_CV`] key, or
    /// when the annotated voltage `cv` satisfies `|cv - target_cv| <
    /// cv_tolerance`. The comparison is strict and has no epsilon, so with the
    /// default tolerance a voltage that differs from the target by exactly
    /// 0.01 in binary arithmetic is dropped. Only the exact key counts; any
    /// other key, such as `FAIMS`, leaves an identification unannotated.
    ///
    /// An integer annotation is converted to `f64` as the source converts an
    /// integer `DataValue` to `double`; integers beyond 2^53 round. A float
    /// annotation is always finite, because
    /// [`MetaValue`](crate::metadata::MetaValue) refuses non-finite floats.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`], without examining any identification,
    /// when `target_cv` is not finite, when `cv_tolerance` is NaN, zero or
    /// negative, and when there are more than
    /// [`MAX_PEPTIDE_IDENTIFICATIONS`](Self::MAX_PEPTIDE_IDENTIFICATIONS)
    /// identifications. The source accepts all of these and silently keeps
    /// only the unannotated identifications, because its strict comparison can
    /// never succeed with them.
    ///
    /// Returns [`Error::InvalidValue`], naming the identification index, when a
    /// `FAIMS_CV` annotation is empty or not a number. The source converts an
    /// empty `DataValue` by throwing `Exception::ConversionError`, and a string
    /// or list `DataValue` by reading an inactive union member (CPP-058).
    ///
    /// Nothing is mutated in any case; on error no partial result escapes.
    pub fn filter_peptides_by_faims_cv(
        peptides: &[PeptideIdentification],
        target_cv: f64,
        cv_tolerance: f64,
    ) -> Result<Vec<PeptideIdentification>> {
        if !target_cv.is_finite() {
            return Err(Error::InvalidValue(
                "the target FAIMS compensation voltage must be finite".into(),
            ));
        }
        if cv_tolerance.is_nan() || cv_tolerance <= 0.0 {
            return Err(Error::InvalidValue(
                "the FAIMS compensation voltage tolerance must be positive".into(),
            ));
        }
        if peptides.len() > Self::MAX_PEPTIDE_IDENTIFICATIONS {
            return Err(Error::InvalidValue(format!(
                "{} peptide identifications exceed the {} this filter allows",
                peptides.len(),
                Self::MAX_PEPTIDE_IDENTIFICATIONS
            )));
        }
        let mut filtered = Vec::new();
        filtered.try_reserve_exact(peptides.len()).map_err(|_| {
            Error::InvalidValue("peptide identification filter allocation failed".into())
        })?;
        for (index, peptide) in peptides.iter().enumerate() {
            let keep = match peptide.metadata.get(FAIMS_CV) {
                None => true,
                Some(value) => {
                    let cv = value.as_f64().map_err(|_| {
                        Error::InvalidValue(format!(
                            "peptide identification {index} has a {FAIMS_CV} annotation that is not a number"
                        ))
                    })?;
                    (cv - target_cv).abs() < cv_tolerance
                }
            };
            if keep {
                filtered.push(peptide.clone());
            }
        }
        Ok(filtered)
    }
}
