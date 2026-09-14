// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Splitting an experiment by FAIMS compensation voltage, the port of
//! `IMDataConverter::splitByFAIMSCV` from `IONMOBILITY/IMDataConverter.h` and
//! `IONMOBILITY/IMDataConverter.cpp`.
//!
//! The source class converts peak maps and spectra between ion-mobility and
//! FAIMS storage models. Only `splitByFAIMSCV` is ported here, the operation
//! FeatureFinderCentroided runs on FAIMS input (decision D5 of the early TOPP
//! bundle). The other members, `reshapeIMFrameToMany`,
//! `splitExperimentByIonMobility`, `reshapeIMFrameToSingle`, `setIMUnit` and
//! `getIMUnit`, are not ported; `splitExperimentByIonMobility` needs the
//! unported `SpectraMerger`.
//!
//! The split is `ImDataConverter::split_by_faims_cv`. It groups the spectra
//! by the voltages `FaimsHelper::get_compensation_voltages` reports, keyed by
//! `FaimsGroupKey`, whose `NotFaims` variant is the source's NaN key made
//! explicit. What the source logs is returned as `FaimsSplitMessage` records,
//! and what it destroys (skipped spectra and, for FAIMS input, the
//! chromatograms) is handed back instead of being dropped silently.
//!
//! Module-level links are plain code spans: the module's `pub mod` line carries
//! its own outer doc comment, and rustdoc resolves combined module docs in the
//! parent module, where these names are not in scope.
//!
//! The work is serial, as the source is. See `docs/IM_DATA_CONVERTER_SUPPORT.md`
//! for the API mapping, the preserved source conventions, the native
//! differences and the evidence.

use super::faims_helper::{CompensationVoltage, FaimsHelper};
use super::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::format::file_info::text_format::{DEFAULT_STREAM_PRECISION, ostream_g};
use crate::metadata::DriftTimeUnit;
use crate::{Error, Result};
use std::fmt;

/// The key of one group of a FAIMS split: the source `double` in each
/// `std::pair<double, MSExperiment>`, with its NaN made explicit.
///
/// The source keys the single group of an experiment without FAIMS
/// compensation voltages with `std::numeric_limits<double>::quiet_NaN()`
/// (`IMDataConverter.cpp:37`). NaN compares unequal to itself and has no place
/// in an ordering, so this port names that case instead:
/// [`NotFaims`](Self::NotFaims). [`volts`](Self::volts) gives back the source
/// value, NaN included.
///
/// The derived order puts `NotFaims` before every voltage and orders the
/// voltages as [`CompensationVoltage`] does. One split never mixes the two
/// variants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaimsGroupKey {
    /// The experiment holds no FAIMS compensation voltage; the group is the
    /// whole, unsplit input. The source key is NaN.
    NotFaims,
    /// The group of one compensation voltage, in volts, with the sign of zero
    /// [`FaimsHelper::get_compensation_voltages`] stored for it.
    Voltage(CompensationVoltage),
}

impl FaimsGroupKey {
    /// The source key: the voltage in volts, or a quiet NaN for
    /// [`NotFaims`](Self::NotFaims).
    pub fn volts(self) -> f64 {
        match self {
            Self::NotFaims => f64::NAN,
            Self::Voltage(voltage) => voltage.volts(),
        }
    }

    /// The voltage, or `None` for [`NotFaims`](Self::NotFaims).
    pub const fn voltage(self) -> Option<CompensationVoltage> {
        match self {
            Self::NotFaims => None,
            Self::Voltage(voltage) => Some(voltage),
        }
    }
}

/// One group of a FAIMS split: the source `std::pair<double, MSExperiment>`.
#[derive(Clone, Debug, PartialEq)]
pub struct FaimsGroup {
    /// The compensation voltage of the group, or
    /// [`FaimsGroupKey::NotFaims`] for the unsplit experiment.
    pub key: FaimsGroupKey,
    /// The spectra of the group in input order.
    ///
    /// For a voltage group the experiment carries a copy of the input's
    /// [`settings`](MSExperiment::settings) and its
    /// [`sql_run_id`](MSExperiment::sql_run_id), and no chromatograms. For
    /// [`FaimsGroupKey::NotFaims`] it is the input experiment itself,
    /// chromatograms included.
    pub experiment: MSExperiment,
}

/// The source log channel a [`FaimsSplitMessage`] is written to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaimsSplitLogLevel {
    /// `OPENMS_LOG_INFO`.
    Info,
    /// `OPENMS_LOG_WARN`.
    Warning,
}

/// A record the source writes to a log stream while splitting, returned
/// instead of logged.
///
/// This kernel module is not wired to a log stream, so the caller decides
/// where the records go; nothing is printed. [`text`](Self::text) (and
/// [`Display`](fmt::Display)) gives the line the source writes.
///
/// There is one record per event. The source log stream collapses repeated
/// identical lines: it prints the first and later reports `<line> occurred N
/// times`. That is presentation of the log stream, not of the split, and it is
/// not reproduced here.
#[derive(Clone, Debug, PartialEq)]
pub enum FaimsSplitMessage {
    /// A warning of [`FaimsHelper::get_compensation_voltages`], verbatim. The
    /// only one it produces is [`FaimsHelper::MISSING_VOLTAGE_WARNING`], when a
    /// FAIMS spectrum carries the unset drift time. Logged first, as the source
    /// logs it from inside `getCompensationVoltages` (`FAIMSHelper.cpp:52`).
    CompensationVoltageWarning(String),
    /// The experiment holds no compensation voltage and is returned unsplit
    /// ([`ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO`], an
    /// `OPENMS_LOG_INFO` record, `IMDataConverter.cpp:35`).
    NoCompensationVoltages,
    /// A FAIMS spectrum whose voltage is not among the detected ones was
    /// skipped (`IMDataConverter.cpp:60`).
    ///
    /// Every FAIMS voltage of the input is detected except the unset sentinel
    /// `-1` ([`ImTypes::DRIFTTIME_NOT_SET`](crate::metadata::ImTypes::DRIFTTIME_NOT_SET)),
    /// which `getCompensationVoltages` removes, so `volts` is always `-1` in
    /// practice.
    UnexpectedCompensationVoltage {
        /// Index of the skipped spectrum in the input.
        spectrum_index: usize,
        /// Its drift time, in volts.
        volts: f64,
    },
    /// A spectrum without a FAIMS voltage was skipped: an MS1 (or MS level 0)
    /// spectrum, or an MS2+ spectrum with no preceding FAIMS spectrum or whose
    /// preceding FAIMS spectrum carried an undetected voltage
    /// (`IMDataConverter.cpp:81`).
    SpectrumWithoutCompensationVoltage {
        /// Index of the skipped spectrum in the input.
        spectrum_index: usize,
    },
}

impl FaimsSplitMessage {
    /// The log channel the source writes the record to.
    pub const fn level(&self) -> FaimsSplitLogLevel {
        match self {
            Self::NoCompensationVoltages => FaimsSplitLogLevel::Info,
            Self::CompensationVoltageWarning(_)
            | Self::UnexpectedCompensationVoltage { .. }
            | Self::SpectrumWithoutCompensationVoltage { .. } => FaimsSplitLogLevel::Warning,
        }
    }

    /// The index of the skipped input spectrum the record reports, if any.
    pub const fn spectrum_index(&self) -> Option<usize> {
        match self {
            Self::UnexpectedCompensationVoltage { spectrum_index, .. }
            | Self::SpectrumWithoutCompensationVoltage { spectrum_index } => Some(*spectrum_index),
            Self::CompensationVoltageWarning(_) | Self::NoCompensationVoltages => None,
        }
    }

    /// The line the source writes, without the trailing newline.
    ///
    /// The voltage of [`UnexpectedCompensationVoltage`](Self::UnexpectedCompensationVoltage)
    /// is formatted as `std::ostream << double` formats it at the log stream's
    /// default precision 6 (`-1` for the sentinel).
    pub fn text(&self) -> String {
        match self {
            Self::CompensationVoltageWarning(text) => text.clone(),
            Self::NoCompensationVoltages => {
                ImDataConverter::NO_COMPENSATION_VOLTAGES_INFO.to_owned()
            }
            Self::UnexpectedCompensationVoltage { volts, .. } => format!(
                "{}{}",
                ImDataConverter::UNEXPECTED_COMPENSATION_VOLTAGE_WARNING,
                ostream_g(*volts, DEFAULT_STREAM_PRECISION)
            ),
            Self::SpectrumWithoutCompensationVoltage { .. } => {
                ImDataConverter::SPECTRUM_WITHOUT_COMPENSATION_VOLTAGE_WARNING.to_owned()
            }
        }
    }
}

impl fmt::Display for FaimsSplitMessage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.text())
    }
}

/// The result of [`ImDataConverter::split_by_faims_cv`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FaimsSplit {
    /// The groups, the source return value.
    ///
    /// Either one group per detected compensation voltage, in ascending voltage
    /// order (`-65`, `-55`, `-45`), or, when the input holds no voltage, exactly
    /// one [`FaimsGroupKey::NotFaims`] group holding the whole input. Never
    /// empty: an empty experiment yields one empty `NotFaims` group, as in the
    /// source.
    pub groups: Vec<FaimsGroup>,
    /// The records the source writes to `OPENMS_LOG_WARN` and
    /// `OPENMS_LOG_INFO`, in the order it writes them.
    pub messages: Vec<FaimsSplitMessage>,
    /// The spectra the source skips, in input order, one for each
    /// [`FaimsSplitMessage::UnexpectedCompensationVoltage`] and
    /// [`FaimsSplitMessage::SpectrumWithoutCompensationVoltage`] record, which
    /// names its input index.
    ///
    /// The source leaves them in the moved-from input and destroys them with
    /// `exp.clear(true)` (`IMDataConverter.cpp:84`); they are returned here so
    /// that nothing is lost silently. Empty for `NotFaims`.
    pub skipped_spectra: Vec<MSSpectrum>,
    /// The chromatograms of a FAIMS input, in input order.
    ///
    /// The source copies only the experimental settings into each voltage group
    /// and destroys the chromatograms with `exp.clear(true)`; no group receives
    /// them. They are returned here instead of being dropped. Empty for
    /// `NotFaims`, whose group keeps them.
    pub dropped_chromatograms: Vec<MSChromatogram>,
}

impl FaimsSplit {
    /// Whether the input held FAIMS compensation voltages, so that the groups
    /// are voltage groups rather than one [`FaimsGroupKey::NotFaims`] group.
    pub fn has_faims(&self) -> bool {
        self.groups
            .iter()
            .any(|group| group.key != FaimsGroupKey::NotFaims)
    }
}

/// Converts peak maps between ion-mobility and FAIMS storage models, the source
/// `class IMDataConverter`.
///
/// The source documentation: *This class converts PeakMaps and MSSpectra
/// from/to different IM/FAIMS storage models.* Only
/// [`split_by_faims_cv`](Self::split_by_faims_cv) is ported.
///
/// The source class holds no state and its members are static; its class test
/// nevertheless constructs and deletes an instance (`IMDataConverter_test.cpp:31-38`).
/// This is therefore a zero-sized unit struct: [`Default`] is the constructor,
/// and it has no drop glue, which is the implicit destructor.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ImDataConverter;

/// Where the split sends one input spectrum.
#[derive(Clone, Copy)]
enum Destination {
    /// Into the voltage group at this index.
    Group(usize),
    /// Skipped: a FAIMS spectrum with an undetected voltage.
    UnexpectedVoltage,
    /// Skipped: no FAIMS voltage and no usable context.
    WithoutVoltage,
}

/// The source's `last_faims_cv` after the latest FAIMS spectrum.
#[derive(Clone, Copy)]
enum Context {
    /// No FAIMS spectrum yet (the source's NaN).
    None,
    /// The latest FAIMS voltage is the group at this index.
    Group(usize),
    /// The latest FAIMS voltage is not a group (the unset sentinel).
    Undetected,
}

/// Reserve exactly `additional` slots, reporting allocation failure as an error
/// instead of aborting.
fn reserve<T>(vector: &mut Vec<T>, additional: usize) -> Result<()> {
    vector
        .try_reserve_exact(additional)
        .map_err(|_| Error::InvalidValue("FAIMS split allocation failed".into()))
}

impl ImDataConverter {
    /// Spectra one [`split_by_faims_cv`](Self::split_by_faims_cv) call may
    /// split.
    ///
    /// Native ceiling; the source splits an experiment of any size. It equals
    /// [`FaimsHelper::MAX_SPECTRA`], the ceiling of the voltage scan the split
    /// starts with, so every experiment whose voltages can be queried can also
    /// be split.
    pub const MAX_SPECTRA: usize = FaimsHelper::MAX_SPECTRA;

    /// The information the source logs when the input holds no compensation
    /// voltage (`IMDataConverter.cpp:35`), verbatim, including its wording.
    pub const NO_COMPENSATION_VOLTAGES_INFO: &'static str =
        "Not FAIMS compensation voltages found in the data. Returning PeakMap as CV NaN.";

    /// The start of the warning the source logs for a FAIMS spectrum whose
    /// voltage was not detected (`IMDataConverter.cpp:60`); the voltage follows.
    pub const UNEXPECTED_COMPENSATION_VOLTAGE_WARNING: &'static str =
        "Encountered spectrum with unexpected FAIMS CV (not in detected set): ";

    /// The warning the source logs for every other skipped spectrum
    /// (`IMDataConverter.cpp:81`), verbatim.
    pub const SPECTRUM_WITHOUT_COMPENSATION_VOLTAGE_WARNING: &'static str =
        "Skipping spectrum without FAIMS CV (no prior FAIMS CV context or unexpected layout).";

    /// Splits an experiment into one experiment per FAIMS compensation voltage
    /// (CV), the source `splitByFAIMSCV(PeakMap&& exp)`.
    ///
    /// The source documentation, carried across: the spectra of the original
    /// experiment are moved to new experiments, so the original is unusable
    /// afterwards.
    ///
    /// - If the dataset contains FAIMS spectra (spectra whose drift time unit is
    ///   [`DriftTimeUnit::FaimsCompensationVoltage`]), the result holds one
    ///   group per unique CV, in ascending CV order. Each group contains all
    ///   spectra assigned to that CV.
    /// - MS2+ spectra without an explicit FAIMS CV are assigned to the last
    ///   seen FAIMS CV, in run order. Spectra without a prior FAIMS CV context
    ///   are skipped with a warning.
    /// - If the dataset contains no FAIMS spectra at all, an informational
    ///   message is logged and a single group is returned. It contains the
    ///   original experiment (no splitting took place) under the key the source
    ///   spells NaN, here [`FaimsGroupKey::NotFaims`]. Callers that
    ///   specifically work on FAIMS data should usually check for FAIMS CVs
    ///   with [`FaimsHelper::get_compensation_voltages`] instead of relying on
    ///   this special case.
    ///
    /// The detailed rules, as the source applies them in one pass over the
    /// spectra in input order (`IMDataConverter.cpp:52-82`):
    ///
    /// - The CVs are those of [`FaimsHelper::get_compensation_voltages`]: every
    ///   FAIMS drift time except the unset sentinel `-1`. Its warnings come
    ///   first in [`FaimsSplit::messages`].
    /// - A FAIMS spectrum, at any MS level, joins the group of its own CV and
    ///   becomes the context for later spectra. If its CV is not a group (only
    ///   the sentinel can be), it is skipped with
    ///   [`FaimsSplitMessage::UnexpectedCompensationVoltage`], and it still
    ///   becomes the context, so the MS2 spectra after it are skipped too.
    /// - A spectrum of any other unit (none, milliseconds, 1/K0, ...) joins the
    ///   context group when its MS level is above 1 and the context is a group.
    ///   Otherwise it is skipped with
    ///   [`FaimsSplitMessage::SpectrumWithoutCompensationVoltage`]. Such a
    ///   spectrum never changes the context: an MS1 spectrum without CV between
    ///   FAIMS spectra is skipped, and the MS2 spectra after it still join the
    ///   earlier CV.
    /// - A CV of `-0.0` and one of `+0.0` are the same group, keyed with the
    ///   sign `get_compensation_voltages` stored; each spectrum keeps its own
    ///   drift time. Infinite CVs are ordinary groups.
    /// - Each voltage group receives a copy of the input's
    ///   [`settings`](MSExperiment::settings) (`IMDataConverter.cpp:46`) and of
    ///   its [`sql_run_id`](MSExperiment::sql_run_id), which the source stores
    ///   as a meta value of the settings. No group receives a chromatogram.
    ///
    /// On success `experiment` is left as [`MSExperiment::default`]: the source
    /// clears the input with `exp.clear(true)` after moving the spectra out,
    /// which also resets its settings, and for input without CVs it moves the
    /// whole experiment into the result. What the source destroys at that point
    /// is returned in [`FaimsSplit::skipped_spectra`] and
    /// [`FaimsSplit::dropped_chromatograms`].
    ///
    /// No ranges are involved: the source builds the groups with `addSpectrum`
    /// and never calls `updateRanges`, so `spectrumRanges().byMSLevel(1)` throws
    /// on every voltage group and FeatureFinderCentroided exits 8 on every FAIMS
    /// input. Ranges here are computed on demand
    /// ([`MSExperiment::spectrum_range_manager`]), so each group's ranges are
    /// always its own. This is a native difference, not emulated.
    ///
    /// Serial and deterministic: groups, spectra and messages follow input
    /// order. Work is one voltage scan, one classification pass and one move
    /// pass; each group's spectrum list is allocated once at its exact size.
    ///
    /// # Errors
    ///
    /// Returns an error, with `experiment` unchanged and nothing moved, when:
    ///
    /// - [`FaimsHelper::get_compensation_voltages`] fails: more than
    ///   [`MAX_SPECTRA`](Self::MAX_SPECTRA) spectra, or a FAIMS spectrum with a
    ///   NaN drift time ([`Error::InvalidValue`] in both cases). The source puts
    ///   the NaN into its `std::set<double>` and `std::map<double, ...>`, whose
    ///   ordering it breaks: executed, a NaN first makes the input unsplit, and a
    ///   later NaN spectrum joins an unrelated group while the MS2 spectra after
    ///   it are skipped;
    /// - the input holds CVs and its settings exceed the resource limits of
    ///   [`ExperimentalSettings::validate`](crate::metadata::ExperimentalSettings::validate),
    ///   checked once before the settings are copied into the groups. The
    ///   source copies settings of any size;
    /// - an allocation fails ([`Error::InvalidValue`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::kernel::im_data_converter::{FaimsGroupKey, FaimsSplitMessage, ImDataConverter};
    /// use openms::metadata::DriftTimeUnit;
    /// use openms::{MSExperiment, MSSpectrum};
    ///
    /// let spectrum = |ms_level: u32, unit: DriftTimeUnit, drift_time: f64| MSSpectrum {
    ///     ms_level,
    ///     drift_time,
    ///     drift_time_unit: unit,
    ///     ..MSSpectrum::default()
    /// };
    /// let faims = DriftTimeUnit::FaimsCompensationVoltage;
    /// let mut experiment = MSExperiment::new();
    /// experiment.spectra = vec![
    ///     spectrum(2, DriftTimeUnit::None, -1.0), // no FAIMS context yet: skipped
    ///     spectrum(1, faims, -45.0),
    ///     spectrum(2, DriftTimeUnit::None, -1.0), // joins -45
    ///     spectrum(1, faims, -65.0),
    /// ];
    ///
    /// let split = ImDataConverter::split_by_faims_cv(&mut experiment)?;
    /// assert!(experiment.is_empty());
    /// let keys: Vec<f64> = split.groups.iter().map(|group| group.key.volts()).collect();
    /// assert_eq!(keys, [-65.0, -45.0]);
    /// assert_eq!(split.groups[1].experiment.len(), 2);
    /// assert_eq!(split.skipped_spectra.len(), 1);
    /// assert_eq!(
    ///     split.messages,
    ///     [FaimsSplitMessage::SpectrumWithoutCompensationVoltage { spectrum_index: 0 }]
    /// );
    ///
    /// let mut plain = MSExperiment::new();
    /// plain.spectra.push(MSSpectrum::default());
    /// let split = ImDataConverter::split_by_faims_cv(&mut plain)?;
    /// assert_eq!(split.groups[0].key, FaimsGroupKey::NotFaims);
    /// assert!(split.groups[0].key.volts().is_nan());
    /// # Ok::<(), openms::Error>(())
    /// ```
    pub fn split_by_faims_cv(experiment: &mut MSExperiment) -> Result<FaimsSplit> {
        let detected = FaimsHelper::get_compensation_voltages(experiment)?;
        let mut messages = Vec::new();
        reserve(&mut messages, detected.warnings.len() + 1)?;
        messages.extend(
            detected
                .warnings
                .into_iter()
                .map(FaimsSplitMessage::CompensationVoltageWarning),
        );

        if detected.voltages.is_empty() {
            messages.push(FaimsSplitMessage::NoCompensationVoltages);
            let mut groups = Vec::new();
            reserve(&mut groups, 1)?;
            groups.push(FaimsGroup {
                key: FaimsGroupKey::NotFaims,
                experiment: std::mem::take(experiment),
            });
            return Ok(FaimsSplit {
                groups,
                messages,
                skipped_spectra: Vec::new(),
                dropped_chromatograms: Vec::new(),
            });
        }

        experiment.settings.validate()?;
        let voltages: Vec<CompensationVoltage> = detected.voltages.into_iter().collect();

        // Classify every spectrum before anything moves, so that a failure
        // leaves the input untouched.
        let mut destinations = Vec::new();
        reserve(&mut destinations, experiment.spectra.len())?;
        let mut counts = Vec::new();
        reserve(&mut counts, voltages.len())?;
        counts.resize(voltages.len(), 0usize);
        let mut skipped = 0usize;
        let mut context = Context::None;
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            let destination = if spectrum.drift_time_unit == DriftTimeUnit::FaimsCompensationVoltage
            {
                // get_compensation_voltages refused NaN, so `new` succeeds; a
                // failure could only mean an undetected voltage.
                let found = CompensationVoltage::new(spectrum.drift_time)
                    .ok()
                    .and_then(|voltage| voltages.binary_search(&voltage).ok());
                match found {
                    Some(group) => {
                        context = Context::Group(group);
                        Destination::Group(group)
                    }
                    None => {
                        context = Context::Undetected;
                        messages.push(FaimsSplitMessage::UnexpectedCompensationVoltage {
                            spectrum_index: index,
                            volts: spectrum.drift_time,
                        });
                        Destination::UnexpectedVoltage
                    }
                }
            } else {
                match context {
                    Context::Group(group) if spectrum.ms_level > 1 => Destination::Group(group),
                    Context::None | Context::Group(_) | Context::Undetected => {
                        messages.push(FaimsSplitMessage::SpectrumWithoutCompensationVoltage {
                            spectrum_index: index,
                        });
                        Destination::WithoutVoltage
                    }
                }
            };
            match destination {
                Destination::Group(group) => counts[group] += 1,
                Destination::UnexpectedVoltage | Destination::WithoutVoltage => skipped += 1,
            }
            destinations.push(destination);
        }

        let mut groups = Vec::new();
        reserve(&mut groups, voltages.len())?;
        for (voltage, &count) in voltages.iter().zip(&counts) {
            let mut spectra = Vec::new();
            reserve(&mut spectra, count)?;
            groups.push(FaimsGroup {
                key: FaimsGroupKey::Voltage(*voltage),
                experiment: MSExperiment {
                    spectra,
                    chromatograms: Vec::new(),
                    settings: experiment.settings.clone(),
                    sql_run_id: experiment.sql_run_id,
                },
            });
        }
        let mut skipped_spectra = Vec::new();
        reserve(&mut skipped_spectra, skipped)?;

        // Commit: nothing below can fail.
        let input = std::mem::take(experiment);
        for (spectrum, destination) in input.spectra.into_iter().zip(destinations) {
            match destination {
                Destination::Group(group) => groups[group].experiment.spectra.push(spectrum),
                Destination::UnexpectedVoltage | Destination::WithoutVoltage => {
                    skipped_spectra.push(spectrum)
                }
            }
        }
        Ok(FaimsSplit {
            groups,
            messages,
            skipped_spectra,
            dropped_chromatograms: input.chromatograms,
        })
    }
}
