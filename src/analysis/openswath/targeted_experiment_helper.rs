// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helper structures shared by the classes of a targeted experiment.
//!
//! Ports `ANALYSIS/TARGETED/TargetedExperimentHelper.h` (with its `.cpp`):
//! the value types that [`ReactionMonitoringTransition`], [`IncludeExcludeTarget`]
//! and [`TargetedExperiment`] are built from, one per TraML element — a
//! controlled vocabulary ([`CV`]), a [`Protein`], a [`Peptide`] or
//! [`Compound`] (both a [`PeptideCompound`]), their [`RetentionTime`]s, a
//! product ion ([`TraMLProduct`]) and its [`Interpretation`]s, plus the
//! bookkeeping records [`Contact`], [`Publication`], [`Instrument`],
//! [`Prediction`] and [`Configuration`]. The two free functions
//! [`get_aa_sequence`] and [`set_modification`] turn a TraML peptide into an
//! [`AASequence`].
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TARGETED_EXPERIMENT_SUPPORT.md)
//! for the member-by-member mapping and the evidence.
//!
//! # Shared conventions
//!
//! * The source derives most records from `CVTermList` or
//!   `CVTermListInterface` (a lazily allocated `CVTermList`). Both become a
//!   `cv_terms: `[`CVTermList`] field here; the lazy allocation is a memory
//!   detail with no observable behaviour, and the `MetaInfoInterface` the
//!   source list carries is [`CVTermList::metadata`].
//! * `Compound` and `Peptide` derive from `PeptideCompound`; here they hold it
//!   as a `base` field and dereference to it, so `peptide.id` and
//!   `peptide.has_retention_time()` read as in the source.
//! * Each source `std::hash` specialisation becomes an [`Hash`] implementation
//!   consistent with [`PartialEq`]. The source hash values themselves are an
//!   implementation detail (an FNV-1a/`hash_combine` mix) that no caller can
//!   rely on, so they are not reproduced; the upstream class test asserts only
//!   that equal values hash equally and that the types work as set keys, which
//!   `tests/targeted_experiment.rs` checks.
//! * Floating-point fields hash with both signs of zero identified, matching
//!   `==`. A value holding NaN is unequal to itself, exactly as in the source,
//!   so it must not be used as a hash-set key; [`Eq`] is implemented on that
//!   documented understanding.
//!
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`ReactionMonitoringTransition`]: crate::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition
//! [`IncludeExcludeTarget`]: crate::analysis::openswath::include_exclude_target::IncludeExcludeTarget
//! [`TargetedExperiment`]: crate::analysis::openswath::targeted_experiment::TargetedExperiment

use crate::chemistry::{AASequence, ModificationsDB, PeptideFragmentType};
use crate::metadata::CVTermList;
use crate::{Error, Result};
use std::hash::{Hash, Hasher};
use std::ops::{Deref, DerefMut};

/// Hash a double so that `0.0` and `-0.0`, which compare equal, hash equally.
pub(crate) fn hash_f64<H: Hasher>(value: f64, state: &mut H) {
    (if value == 0.0 { 0 } else { value.to_bits() }).hash(state);
}

/// An instrument configuration used to acquire a product ion or a target.
///
/// TraML `<Configuration>`: CV terms describing the configuration, references
/// to the contact and the instrument, and a list of validation records.
///
/// Source `Configuration` declares no `operator==`, so the inherited
/// `CVTermList::operator==` compares only the CV terms and ignores
/// `contact_ref`, `instrument_ref` and `validations`. This port compares every
/// field: two configurations naming different instruments are different. The
/// source behaviour leaks into `TraMLProduct` and `IncludeExcludeTarget`
/// equality; it is recorded as a C++ defect candidate in the support document.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Configuration {
    /// CV terms of the configuration (the source base class).
    pub cv_terms: CVTermList,
    /// Reference to a [`Contact`] by id.
    pub contact_ref: String,
    /// Reference to an [`Instrument`] by id.
    pub instrument_ref: String,
    /// Validation records, each a CV term list.
    pub validations: Vec<CVTermList>,
}

impl Eq for Configuration {}

/// A controlled vocabulary referenced by a TraML document (`<cv>`).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct CV {
    /// Short identifier, as used in the `cvRef` of a CV term (e.g. `MS`).
    pub id: String,
    /// Full name of the vocabulary.
    pub fullname: String,
    /// Version of the vocabulary.
    pub version: String,
    /// Location of the vocabulary. The source member is `URI`.
    pub uri: String,
}

impl CV {
    /// A vocabulary record from its four fields, as the source's only
    /// constructor. The source has no default constructor; [`Default`] is a
    /// native addition giving four empty strings.
    pub fn new(
        id: impl Into<String>,
        fullname: impl Into<String>,
        version: impl Into<String>,
        uri: impl Into<String>,
    ) -> Self {
        Self {
            id: id.into(),
            fullname: fullname.into(),
            version: version.into(),
            uri: uri.into(),
        }
    }
}

/// A protein referenced by the peptides of a targeted experiment.
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Protein {
    /// CV terms of the protein (the source base class).
    pub cv_terms: CVTermList,
    /// Unique identifier; peptides refer to it through
    /// [`Peptide::protein_refs`].
    pub id: String,
    /// Amino acid sequence, if known.
    pub sequence: String,
}

impl Eq for Protein {}

/// Unit of a stored [`RetentionTime`].
///
/// Discriminants are the source's `std::int8_t` values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(i8)]
pub enum RTUnit {
    /// Retention time stored in seconds.
    Second = 0,
    /// Retention time stored in minutes.
    Minute = 1,
    /// No stored annotation.
    Unknown = 2,
    /// The source's `SIZE_OF_RTUNIT` enumerator. It is also what a
    /// default-constructed `RetentionTime` holds, so readers of the source use
    /// it to mean "no unit was assigned at all", as opposed to [`Self::Unknown`].
    SizeOfRTUnit = 3,
}

/// Kind of a stored [`RetentionTime`].
///
/// Discriminants are the source's `std::int8_t` values.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(i8)]
pub enum RTType {
    /// Undefined local chromatography.
    Local = 0,
    /// Standardised reference chromatography.
    Normalized = 1,
    /// Predicted by referenced software.
    Predicted = 2,
    /// H-PINS, "the de facto standard providing the retention times".
    Hpins = 3,
    /// iRT retention time standard.
    Irt = 4,
    /// No stored annotation.
    Unknown = 5,
    /// The source's `SIZE_OF_RTTYPE` enumerator, which a default-constructed
    /// `RetentionTime` holds: no type was assigned at all.
    SizeOfRTType = 6,
}

/// A retention time, as stored in a TraML `<RetentionTime>` element.
///
/// According to the standard each retention time can carry one or more CV
/// terms describing it. The unit and type are stored in
/// [`retention_time_unit`](Self::retention_time_unit) and
/// [`retention_time_type`](Self::retention_time_type); the value itself is
/// read with [`rt`](Self::rt) and written with [`set_rt`](Self::set_rt).
/// Retention-time windows and lower/upper limits are not modelled as fields
/// (the source leaves them commented out) but remain expressible as CV terms.
#[derive(Clone, Debug, PartialEq)]
pub struct RetentionTime {
    /// CV terms of the retention time (the source `CVTermListInterface` base).
    pub cv_terms: CVTermList,
    /// Reference to the software that produced the value.
    pub software_ref: String,
    /// Unit; defaults to [`RTUnit::SizeOfRTUnit`], as the source.
    pub retention_time_unit: RTUnit,
    /// Kind; defaults to [`RTType::SizeOfRTType`], as the source.
    pub retention_time_type: RTType,
    retention_time_set: bool,
    retention_time: f64,
}

impl Default for RetentionTime {
    fn default() -> Self {
        Self {
            cv_terms: CVTermList::default(),
            software_ref: String::new(),
            retention_time_unit: RTUnit::SizeOfRTUnit,
            retention_time_type: RTType::SizeOfRTType,
            retention_time_set: false,
            retention_time: 0.0,
        }
    }
}

impl Eq for RetentionTime {}

impl Hash for RetentionTime {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cv_terms.hash(state);
        self.software_ref.hash(state);
        self.retention_time_unit.hash(state);
        self.retention_time_type.hash(state);
        self.retention_time_set.hash(state);
        hash_f64(self.retention_time, state);
    }
}

impl RetentionTime {
    /// A retention time with no value, unit or type assigned.
    pub fn new() -> Self {
        Self::default()
    }

    /// Whether a value has been stored with [`set_rt`](Self::set_rt).
    pub fn is_rt_set(&self) -> bool {
        self.retention_time_set
    }

    /// Store the retention time value and mark it as set.
    pub fn set_rt(&mut self, rt: f64) {
        self.retention_time = rt;
        self.retention_time_set = true;
    }

    /// The stored retention time, or `None` when none was set.
    ///
    /// Source `getRT` states `isRTset()` as an `OPENMS_PRECONDITION`, which
    /// release builds do not check (they return `0.0`); the unset case is
    /// `None` here.
    pub fn rt(&self) -> Option<f64> {
        self.retention_time_set.then_some(self.retention_time)
    }
}

/// Base of [`Peptide`] and [`Compound`]: identifier, retention times, charge
/// and precursor ion-mobility drift time.
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideCompound {
    /// CV terms of the analyte (the source base class).
    pub cv_terms: CVTermList,
    /// Unique identifier; transitions refer to it by peptide or compound
    /// reference.
    pub id: String,
    /// Retention times; only the first is read by the retention-time
    /// accessors below.
    pub rts: Vec<RetentionTime>,
    charge: Option<i32>,
    drift_time: f64,
}

impl Default for PeptideCompound {
    fn default() -> Self {
        Self {
            cv_terms: CVTermList::default(),
            id: String::new(),
            rts: Vec::new(),
            charge: None,
            drift_time: -1.0,
        }
    }
}

/// Source `operator==` compares the CV terms, `rts`, `id` and the charge but
/// **not** the drift time; this port keeps that, so two analytes differing
/// only in drift time compare (and hash) equal.
impl Eq for PeptideCompound {}

impl Hash for PeptideCompound {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cv_terms.hash(state);
        self.rts.hash(state);
        self.id.hash(state);
        self.charge.hash(state);
    }
}

impl PeptideCompound {
    /// An analyte with no id, retention time or charge, and drift time `-1`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the charge state and mark it as set.
    pub fn set_charge_state(&mut self, charge: i32) {
        self.charge = Some(charge);
    }

    /// Whether a charge state has been set.
    pub fn has_charge(&self) -> bool {
        self.charge.is_some()
    }

    /// The charge state, or `None` when none was set.
    ///
    /// Source `getChargeState` guards this with an `OPENMS_PRECONDITION` that
    /// release builds skip, returning `0`; the unset case is `None` here.
    pub fn charge_state(&self) -> Option<i32> {
        self.charge
    }

    /// Set the precursor ion-mobility drift time. `-1` is the source's "unset"
    /// default; no value is rejected.
    pub fn set_drift_time(&mut self, dt: f64) {
        self.drift_time = dt;
    }

    /// The precursor ion-mobility drift time; `-1` when never set.
    pub fn drift_time(&self) -> f64 {
        self.drift_time
    }

    /// Whether the analyte has an annotated retention time: the first entry
    /// of [`rts`](Self::rts) exists and has a value set.
    pub fn has_retention_time(&self) -> bool {
        self.rts.first().is_some_and(RetentionTime::is_rt_set)
    }

    fn first_rt(&self) -> Result<&RetentionTime> {
        match self.rts.first() {
            Some(rt) if rt.is_rt_set() => Ok(rt),
            _ => Err(Error::InvalidValue(
                "No retention time information available".into(),
            )),
        }
    }

    /// The retention time value of the first entry of [`rts`](Self::rts).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] (the source's `Exception::IllegalArgument`,
    /// with its message) when [`has_retention_time`](Self::has_retention_time)
    /// is false.
    pub fn retention_time(&self) -> Result<f64> {
        let rt = self.first_rt()?;
        Ok(rt.retention_time)
    }

    /// The retention time type of the first entry of [`rts`](Self::rts).
    ///
    /// # Errors
    ///
    /// As [`retention_time`](Self::retention_time).
    pub fn retention_time_type(&self) -> Result<RTType> {
        Ok(self.first_rt()?.retention_time_type)
    }

    /// The retention time unit (minutes or seconds) of the first entry of
    /// [`rts`](Self::rts).
    ///
    /// # Errors
    ///
    /// As [`retention_time`](Self::retention_time).
    pub fn retention_time_unit(&self) -> Result<RTUnit> {
        Ok(self.first_rt()?.retention_time_unit)
    }
}

/// A small molecule: a [`PeptideCompound`] with a theoretical mass, a SMILES
/// string and a molecular formula.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Compound {
    /// Identifier, retention times, charge and drift time.
    pub base: PeptideCompound,
    /// Molecular (sum) formula, as text.
    pub molecular_formula: String,
    /// SMILES string.
    pub smiles_string: String,
    /// Theoretical (neutral) mass; `0` by default.
    pub theoretical_mass: f64,
}

impl Eq for Compound {}

impl Hash for Compound {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.base.hash(state);
        self.molecular_formula.hash(state);
        self.smiles_string.hash(state);
        hash_f64(self.theoretical_mass, state);
    }
}

impl Deref for Compound {
    type Target = PeptideCompound;
    fn deref(&self) -> &PeptideCompound {
        &self.base
    }
}

impl DerefMut for Compound {
    fn deref_mut(&mut self) -> &mut PeptideCompound {
        &mut self.base
    }
}

/// A modification on a TraML [`Peptide`] (source `Peptide::Modification`).
///
/// `location` is a residue index, with `-1` for the N-terminus and the
/// sequence length for the C-terminus (see [`set_modification`]).
///
/// The source default constructor leaves `avg_mass_delta` and
/// `mono_mass_delta` uninitialised; they are `0` here. The source declares no
/// `operator==`, so the inherited `CVTermListInterface::operator==` compares
/// only the CV terms and ignores the masses, location and UniMod id; this port
/// compares every field (C++ defect candidate, see the support document).
#[derive(Clone, Debug, PartialEq)]
pub struct PeptideModification {
    /// CV terms of the modification.
    pub cv_terms: CVTermList,
    /// Average mass difference in Da.
    pub avg_mass_delta: f64,
    /// Monoisotopic mass difference in Da.
    pub mono_mass_delta: f64,
    /// Residue index; `-1` N-terminus, sequence length C-terminus.
    pub location: i32,
    /// UniMod record id, `-1` when unknown.
    pub unimod_id: i32,
}

impl Default for PeptideModification {
    fn default() -> Self {
        Self {
            cv_terms: CVTermList::default(),
            avg_mass_delta: 0.0,
            mono_mass_delta: 0.0,
            location: -1,
            unimod_id: -1,
        }
    }
}

impl Eq for PeptideModification {}

impl Hash for PeptideModification {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cv_terms.hash(state);
        hash_f64(self.avg_mass_delta, state);
        hash_f64(self.mono_mass_delta, state);
        self.location.hash(state);
        self.unimod_id.hash(state);
    }
}

/// A peptide: a [`PeptideCompound`] with its unmodified sequence, its
/// modifications and the proteins it belongs to.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Peptide {
    /// Identifier, retention times, charge and drift time.
    pub base: PeptideCompound,
    /// Ids of the [`Protein`]s this peptide maps to.
    pub protein_refs: Vec<String>,
    /// Evidence for the peptide, as CV terms.
    pub evidence: CVTermList,
    /// The "naked" amino acid sequence without modifications.
    pub sequence: String,
    /// Modifications on [`sequence`](Self::sequence).
    pub mods: Vec<PeptideModification>,
    peptide_group_label: String,
}

impl Eq for Peptide {}

impl Hash for Peptide {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.base.hash(state);
        self.protein_refs.hash(state);
        self.evidence.hash(state);
        self.sequence.hash(state);
        self.mods.hash(state);
        self.peptide_group_label.hash(state);
    }
}

impl Deref for Peptide {
    type Target = PeptideCompound;
    fn deref(&self) -> &PeptideCompound {
        &self.base
    }
}

impl DerefMut for Peptide {
    fn deref_mut(&mut self) -> &mut PeptideCompound {
        &mut self.base
    }
}

impl Peptide {
    /// Set the peptide group label.
    ///
    /// The label names the non-labelled peptide group the peptide belongs to
    /// (PSI-MS `MS:1000893`): an arbitrary string marking a set of peptides
    /// that belong together and differ only by isotopic labels, so the heavy
    /// and light forms of one peptide share a label.
    pub fn set_peptide_group_label(&mut self, label: impl Into<String>) {
        self.peptide_group_label = label.into();
    }

    /// The peptide group label; empty when unset.
    pub fn peptide_group_label(&self) -> &str {
        &self.peptide_group_label
    }
}

/// A contact person or organisation (`<Contact>`).
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Contact {
    /// CV terms describing the contact.
    pub cv_terms: CVTermList,
    /// Unique identifier.
    pub id: String,
}

impl Eq for Contact {}

/// A publication (`<Publication>`).
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Publication {
    /// CV terms describing the publication.
    pub cv_terms: CVTermList,
    /// Unique identifier.
    pub id: String,
}

impl Eq for Publication {}

/// An instrument (`<Instrument>`).
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Instrument {
    /// CV terms describing the instrument.
    pub cv_terms: CVTermList,
    /// Unique identifier.
    pub id: String,
}

impl Eq for Instrument {}

/// A prediction of a suitable transition by some software (`<Prediction>`).
#[derive(Clone, Debug, Default, PartialEq, Hash)]
pub struct Prediction {
    /// CV terms of the prediction.
    pub cv_terms: CVTermList,
    /// Reference to the predicting software.
    pub software_ref: String,
    /// Reference to the responsible contact.
    pub contact_ref: String,
}

impl Eq for Prediction {}

/// Interpretation of a product ion: ion type, series ordinal and rank.
///
/// The source keeps a commented-out ion-type list with PSI-MS accessions; the
/// type used is `Residue::ResidueType`, which is [`PeptideFragmentType`] here.
#[derive(Clone, Debug, PartialEq, Hash)]
pub struct Interpretation {
    /// CV terms of the interpretation.
    pub cv_terms: CVTermList,
    /// `MS:1000903`, product ion series ordinal (e.g. 8 for a y8 ion).
    pub ordinal: u8,
    /// `MS:1000926`, product interpretation rank (e.g. 1 for the most likely).
    pub rank: u8,
    /// Ion type (b, y, z, ...). Defaults to
    /// [`PeptideFragmentType::Unannotated`], which implies no PSI-MS term.
    pub iontype: PeptideFragmentType,
}

impl Eq for Interpretation {}

impl Default for Interpretation {
    fn default() -> Self {
        Self {
            cv_terms: CVTermList::default(),
            ordinal: 0,
            rank: 0,
            iontype: PeptideFragmentType::Unannotated,
        }
    }
}

/// A product ion (`<Product>` or `<IntermediateProduct>` in TraML).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TraMLProduct {
    /// CV terms of the product ion.
    pub cv_terms: CVTermList,
    charge: Option<i32>,
    mz: f64,
    configuration_list: Vec<Configuration>,
    interpretation_list: Vec<Interpretation>,
}

impl Eq for TraMLProduct {}

impl Hash for TraMLProduct {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.cv_terms.hash(state);
        self.charge.hash(state);
        hash_f64(self.mz, state);
        self.configuration_list.hash(state);
        self.interpretation_list.hash(state);
    }
}

impl TraMLProduct {
    /// A product with m/z `0`, no charge and no configuration or
    /// interpretation.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the product ion charge and mark it as set.
    pub fn set_charge_state(&mut self, charge: i32) {
        self.charge = Some(charge);
    }

    /// Whether the product charge has been set.
    pub fn has_charge(&self) -> bool {
        self.charge.is_some()
    }

    /// The product charge, or `None` when none was set (the source's
    /// unchecked precondition returns `0`).
    pub fn charge_state(&self) -> Option<i32> {
        self.charge
    }

    /// The product ion m/z; `0` by default.
    pub fn mz(&self) -> f64 {
        self.mz
    }

    /// Set the product ion m/z. No value is rejected, as in the source.
    pub fn set_mz(&mut self, mz: f64) {
        self.mz = mz;
    }

    /// The instrument configurations used for this product.
    pub fn configuration_list(&self) -> &[Configuration] {
        &self.configuration_list
    }

    /// Append an instrument configuration.
    pub fn add_configuration(&mut self, configuration: Configuration) {
        self.configuration_list.push(configuration);
    }

    /// The interpretations of this product ion.
    pub fn interpretation_list(&self) -> &[Interpretation] {
        &self.interpretation_list
    }

    /// Append an interpretation.
    pub fn add_interpretation(&mut self, interpretation: Interpretation) {
        self.interpretation_list.push(interpretation);
    }

    /// Remove all interpretations.
    pub fn reset_interpretations(&mut self) {
        self.interpretation_list.clear();
    }
}

/// Set `modification` on `aas` at a TraML location.
///
/// `location == -1` sets the N-terminal modification, `location == max_size`
/// the C-terminal one, any other value the residue at that index. `max_size`
/// is the sequence length.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `location` is outside `-1..=max_size` (the
/// source states this as an `OPENMS_PRECONDITION` that release builds skip,
/// after which `AASequence::setModification` receives a wrapped index), and
/// whatever [`AASequence::set_modification`] and its terminal counterparts
/// return for an unknown modification name or a residue index past the end
/// of `aas`. On error `aas` is unchanged.
pub fn set_modification(
    location: i32,
    max_size: i32,
    modification: &str,
    aas: &mut AASequence,
) -> Result<()> {
    set_modification_with_registry(
        location,
        max_size,
        modification,
        aas,
        ModificationsDB::global(),
    )
}

/// [`set_modification`] against a caller-owned modification registry.
///
/// # Errors
///
/// As [`set_modification`].
pub fn set_modification_with_registry(
    location: i32,
    max_size: i32,
    modification: &str,
    aas: &mut AASequence,
    db: &ModificationsDB,
) -> Result<()> {
    if location < -1 || location > max_size {
        return Err(Error::InvalidValue(format!(
            "Location has invalid value {location}"
        )));
    }
    if location == -1 {
        aas.set_n_terminal_modification_with_registry(modification, db)
    } else if location == max_size {
        aas.set_c_terminal_modification_with_registry(modification, db)
    } else {
        let index = usize::try_from(location)
            .map_err(|_| Error::InvalidValue(format!("Location has invalid value {location}")))?;
        aas.set_modification_with_registry(index, modification, db)
    }
}

/// What [`get_aa_sequence_with`] does with a modification it cannot resolve.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum UnresolvedModificationPolicy {
    /// Refuse the conversion with [`Error::InvalidValue`]. The native default,
    /// because dropping a modification silently changes the analyte.
    #[default]
    Reject,
    /// Reproduce the source: skip the modification, keep converting, and
    /// report it as a warning (the source prints it to `stderr`).
    SourceSkip,
}

/// Result of [`get_aa_sequence_with`]: the sequence and the warnings the
/// source would have logged.
#[derive(Clone, Debug, PartialEq)]
pub struct AASequenceConversion {
    /// The modified sequence.
    pub sequence: AASequence,
    /// One line per warning, in the order the source logs them: a missing
    /// UniMod id (which triggers the mass lookup) and, under
    /// [`UnresolvedModificationPolicy::SourceSkip`], each skipped
    /// modification.
    pub warnings: Vec<String>,
}

/// Convert a TraML [`Peptide`] to an [`AASequence`] with its modifications.
///
/// Uses the global [`ModificationsDB`] and
/// [`UnresolvedModificationPolicy::Reject`]; warnings are discarded. See
/// [`get_aa_sequence_with`] for the procedure.
///
/// # Errors
///
/// As [`get_aa_sequence_with`] under `Reject`.
pub fn get_aa_sequence(peptide: &Peptide) -> Result<AASequence> {
    get_aa_sequence_with(
        peptide,
        ModificationsDB::global(),
        UnresolvedModificationPolicy::Reject,
    )
    .map(|conversion| conversion.sequence)
}

/// Convert a TraML [`Peptide`] to an [`AASequence`], choosing the registry
/// and the handling of unresolvable modifications.
///
/// [`Peptide::sequence`] is the naked sequence; it is parsed first, then each
/// entry of [`Peptide::mods`] is applied in order:
///
/// 1. When `unimod_id != -1`, the modification `UniMod:<id>` is set at
///    `location` with [`set_modification`].
/// 2. Otherwise a warning is recorded and the modification is looked up by
///    its monoisotopic mass difference within 1.0 Da on the residue at
///    `location` (the closest match strictly inside the tolerance, as the
///    source `getBestModificationByDiffMonoMass`), and its id is set.
/// 3. When no modification matches, the policy decides.
///
/// For a terminal location (`-1` or the sequence length) without a UniMod id
/// there is no residue to match: the source then reads
/// `peptide.sequence[location]`, which is undefined behaviour for `-1` and
/// the terminating NUL for the length. This port treats such a modification
/// as unresolvable rather than guessing a residue.
///
/// # Errors
///
/// [`Error::InvalidValue`] or [`Error::Parse`] when the sequence does not
/// parse, when a location is outside `-1..=len`, when a UniMod id is unknown
/// to `db`, and, under [`UnresolvedModificationPolicy::Reject`], when a
/// modification cannot be resolved by mass. The source throws for the first
/// three as well (from `AASequence`) and skips the last.
pub fn get_aa_sequence_with(
    peptide: &Peptide,
    db: &ModificationsDB,
    policy: UnresolvedModificationPolicy,
) -> Result<AASequenceConversion> {
    let mut aas = AASequence::parse_with_registry(&peptide.sequence, db)?;
    let max_size = i32::try_from(peptide.sequence.len()).map_err(|_| {
        Error::InvalidValue("peptide sequence is too long for a TraML location".into())
    })?;
    let mut warnings = Vec::new();
    for modification in &peptide.mods {
        if modification.unimod_id != -1 {
            let name = format!("UniMod:{}", modification.unimod_id);
            set_modification_with_registry(modification.location, max_size, &name, &mut aas, db)?;
            continue;
        }
        warnings.push(format!(
            "Warning: No UniMod id set for modification on peptide {}. Will try to infer modification id by mass next.",
            peptide.sequence
        ));
        let residue = usize::try_from(modification.location)
            .ok()
            .and_then(|index| peptide.sequence.as_bytes().get(index))
            .map(|&byte| char::from(byte));
        let found = match residue {
            Some(residue) => {
                db.best_by_mass(modification.mono_mass_delta, 1.0, Some(residue), None)?
            }
            None => None,
        };
        match found {
            Some(found) => {
                let id = found.name().to_owned();
                set_modification_with_registry(modification.location, max_size, &id, &mut aas, db)?;
            }
            None => {
                let message = format!(
                    "Warning: Could not determine modification with delta mass {} for peptide {} at position {}",
                    modification.mono_mass_delta, peptide.sequence, modification.location
                );
                match policy {
                    UnresolvedModificationPolicy::Reject => {
                        return Err(Error::InvalidValue(message));
                    }
                    UnresolvedModificationPolicy::SourceSkip => {
                        warnings.push(message);
                        warnings.push("Skipping this modification".into());
                    }
                }
            }
        }
    }
    Ok(AASequenceConversion {
        sequence: aas,
        warnings,
    })
}
