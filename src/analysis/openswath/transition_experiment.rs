// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The light (memory-lean) targeted experiment OpenSWATH works on: transitions,
//! compounds (peptides or metabolites) and proteins.
//!
//! Ports `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h` (header-only),
//! including its `std::hash` specialisations. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md)
//! for the member mapping and evidence.
//!
//! The source's typedef members `LightTargetedExperiment::Transition`,
//! `::Peptide`, `::Compound` and `::Protein` have no Rust counterpart
//! (inherent associated types are unstable); they name [`LightTransition`](crate::analysis::openswath::transition_experiment::LightTransition),
//! [`LightCompound`](crate::analysis::openswath::transition_experiment::LightCompound), [`LightCompound`](crate::analysis::openswath::transition_experiment::LightCompound) and [`LightProtein`](crate::analysis::openswath::transition_experiment::LightProtein).
//!
//! # Identity
//!
//! As in the source, equality and hashing of [`LightTransition`](crate::analysis::openswath::transition_experiment::LightTransition),
//! [`LightCompound`](crate::analysis::openswath::transition_experiment::LightCompound) and [`LightProtein`](crate::analysis::openswath::transition_experiment::LightProtein) look at the identifier only (the
//! transition name, the compound id, the protein id), so two records with the
//! same identifier and different payloads compare equal. The hash values
//! differ from the source's `std::hash` results; they are process-local in
//! both and never persisted.

use std::collections::HashMap;
use std::fmt;
use std::hash::{Hash, Hasher};

/// Compact enum for fragment ion types, replacing string storage (one byte
/// instead of a `std::string`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum FragmentIonType {
    /// A type string that none of the other variants names.
    Unknown = 0,
    /// a-ion.
    AIon,
    /// b-ion.
    BIon,
    /// c-ion.
    CIon,
    /// x-ion.
    XIon,
    /// y-ion.
    YIon,
    /// z-ion.
    ZIon,
    /// z'-ion (z prime).
    ZPrimeIon,
    /// z.-ion (z dot).
    ZDotIon,
    /// Precursor ion.
    Precursor,
    /// b-ion with water loss.
    BMinusH2O,
    /// y-ion with water loss.
    YMinusH2O,
    /// b-ion with ammonia loss.
    BMinusNH3,
    /// y-ion with ammonia loss.
    YMinusNH3,
    /// No fragment type set (source value 255).
    #[default]
    Empty = 255,
}

impl FragmentIonType {
    /// Convert a fragment ion type string (e.g. `"y"`, `"b"`, `"prec"`) to the
    /// enum. Source `stringToFragmentIonType`.
    ///
    /// The empty string is [`Self::Empty`]; `"prec"` and `"precursor"` are
    /// [`Self::Precursor`]; the misspellings `"b-H20"` and `"y-H20"` (digit
    /// zero) are accepted as water losses, as in the source; anything
    /// unrecognised is [`Self::Unknown`]. Matching is case-sensitive.
    pub fn from_type_str(s: &str) -> Self {
        match s {
            "" => Self::Empty,
            "a" => Self::AIon,
            "b" => Self::BIon,
            "c" => Self::CIon,
            "x" => Self::XIon,
            "y" => Self::YIon,
            "z" => Self::ZIon,
            "z'" => Self::ZPrimeIon,
            "z." => Self::ZDotIon,
            "prec" | "precursor" => Self::Precursor,
            "b-H2O" | "b-H20" => Self::BMinusH2O,
            "y-H2O" | "y-H20" => Self::YMinusH2O,
            "b-NH3" => Self::BMinusNH3,
            "y-NH3" => Self::YMinusNH3,
            _ => Self::Unknown,
        }
    }

    /// The fragment type string (e.g. `"y"`, `"b"`, `"prec"`). Source
    /// `fragmentIonTypeToString`.
    ///
    /// [`Self::Empty`] is `""` and [`Self::Unknown`] is `"unknown"`, so
    /// `from_type_str(as_str())` round-trips every variant except `Unknown`
    /// (which reads back as `Unknown` only because `"unknown"` is itself
    /// unrecognised).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Empty => "",
            Self::AIon => "a",
            Self::BIon => "b",
            Self::CIon => "c",
            Self::XIon => "x",
            Self::YIon => "y",
            Self::ZIon => "z",
            Self::ZPrimeIon => "z'",
            Self::ZDotIon => "z.",
            Self::Precursor => "prec",
            Self::BMinusH2O => "b-H2O",
            Self::YMinusH2O => "y-H2O",
            Self::BMinusNH3 => "b-NH3",
            Self::YMinusNH3 => "y-NH3",
            Self::Unknown => "unknown",
        }
    }
}

impl fmt::Display for FragmentIonType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Packed boolean flags of a transition, one byte instead of four `bool`s.
///
/// Defaults, as in the source: not a decoy, detecting, quantifying, not
/// identifying. The source's four reserved bits are kept zero.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct TransitionFlags(u8);

impl TransitionFlags {
    const DECOY: u8 = 1;
    const DETECTING: u8 = 1 << 1;
    const QUANTIFYING: u8 = 1 << 2;
    const IDENTIFYING: u8 = 1 << 3;

    fn get(self, bit: u8) -> bool {
        self.0 & bit != 0
    }

    fn set(&mut self, bit: u8, value: bool) {
        if value {
            self.0 |= bit;
        } else {
            self.0 &= !bit;
        }
    }

    /// Whether the transition is a decoy.
    pub fn decoy(self) -> bool {
        self.get(Self::DECOY)
    }
    /// Set the decoy flag.
    pub fn set_decoy(&mut self, value: bool) {
        self.set(Self::DECOY, value);
    }
    /// Whether the transition is used for detection.
    pub fn detecting(self) -> bool {
        self.get(Self::DETECTING)
    }
    /// Set the detecting flag.
    pub fn set_detecting(&mut self, value: bool) {
        self.set(Self::DETECTING, value);
    }
    /// Whether the transition is used for quantification.
    pub fn quantifying(self) -> bool {
        self.get(Self::QUANTIFYING)
    }
    /// Set the quantifying flag.
    pub fn set_quantifying(&mut self, value: bool) {
        self.set(Self::QUANTIFYING, value);
    }
    /// Whether the transition is used for identification (IPF).
    pub fn identifying(self) -> bool {
        self.get(Self::IDENTIFYING)
    }
    /// Set the identifying flag.
    pub fn set_identifying(&mut self, value: bool) {
        self.set(Self::IDENTIFYING, value);
    }
}

impl Default for TransitionFlags {
    fn default() -> Self {
        Self(Self::DETECTING | Self::QUANTIFYING)
    }
}

impl fmt::Debug for TransitionFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("TransitionFlags")
            .field("decoy", &self.decoy())
            .field("detecting", &self.detecting())
            .field("quantifying", &self.quantifying())
            .field("identifying", &self.identifying())
            .finish()
    }
}

/// A light transition: one precursor/product m/z pair of an assay.
///
/// Equality and hashing use [`Self::transition_name`] only (see the module notes).
#[derive(Clone, Debug)]
pub struct LightTransition {
    /// Transition name, the native ID of its chromatogram.
    pub transition_name: String,
    /// Id of the compound (peptide or metabolite) this transition belongs to.
    pub peptide_ref: String,
    /// Library (expected relative) intensity.
    pub library_intensity: f64,
    /// Product (fragment) m/z.
    pub product_mz: f64,
    /// Precursor m/z.
    pub precursor_mz: f64,
    /// Precursor ion mobility; -1 means unset.
    pub precursor_im: f64,
    /// Fragment charge (compact: the range is typically 1-8); 0 means unset.
    pub fragment_charge: i8,
    /// Packed boolean flags.
    pub flags: TransitionFlags,
    /// Fragment ion ordinal (e.g. 7 for y7), for round-trip TSV/PQP I/O; -1
    /// means unset.
    pub fragment_nr: i16,
    /// Fragment ion type, for round-trip TSV/PQP I/O.
    pub fragment_type: FragmentIonType,
    /// Peptidoforms, for IPF.
    pub peptidoforms: Vec<String>,
}

impl Default for LightTransition {
    fn default() -> Self {
        Self {
            transition_name: String::new(),
            peptide_ref: String::new(),
            library_intensity: 0.0,
            product_mz: 0.0,
            precursor_mz: 0.0,
            precursor_im: -1.0,
            fragment_charge: 0,
            flags: TransitionFlags::default(),
            fragment_nr: -1,
            fragment_type: FragmentIonType::Empty,
            peptidoforms: Vec::new(),
        }
    }
}

impl PartialEq for LightTransition {
    fn eq(&self, other: &Self) -> bool {
        self.transition_name == other.transition_name
    }
}

impl Eq for LightTransition {}

impl Hash for LightTransition {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.transition_name.hash(state);
    }
}

impl LightTransition {
    /// Whether the transition is a decoy (legacy accessor over [`Self::flags`]).
    pub fn decoy(&self) -> bool {
        self.flags.decoy()
    }

    /// Set the decoy flag.
    pub fn set_decoy(&mut self, decoy: bool) {
        self.flags.set_decoy(decoy);
    }

    /// The fragment type as a string, for backward compatibility; see
    /// [`FragmentIonType::as_str`].
    pub fn fragment_type_str(&self) -> &'static str {
        self.fragment_type.as_str()
    }

    /// Set the fragment type from a string; see [`FragmentIonType::from_type_str`].
    pub fn set_fragment_type(&mut self, s: &str) {
        self.fragment_type = FragmentIonType::from_type_str(s);
    }

    /// Reconstruct the annotation from fragment type, number and charge.
    ///
    /// Returns strings like `"b8"`, `"y6^2"` or `"prec"`, or the empty string
    /// when no fragment type is set. The number is appended when it is
    /// non-negative, the charge (as `^z`) when it is positive; an unknown type
    /// yields `"unknown"` followed by those parts.
    pub fn annotation(&self) -> String {
        match self.fragment_type {
            FragmentIonType::Empty => String::new(),
            FragmentIonType::Precursor => "prec".to_string(),
            other => {
                let mut result = other.as_str().to_string();
                if self.fragment_nr >= 0 {
                    result.push_str(&self.fragment_nr.to_string());
                }
                if self.fragment_charge > 0 {
                    result.push('^');
                    result.push_str(&self.fragment_charge.to_string());
                }
                result
            }
        }
    }

    /// The product (fragment) charge state; 0 when unset.
    pub fn product_charge_state(&self) -> i32 {
        i32::from(self.fragment_charge)
    }

    /// Whether a product charge is set (non-zero).
    pub fn is_product_charge_state_set(&self) -> bool {
        self.fragment_charge != 0
    }

    /// Whether a precursor ion mobility is set, i.e. it is not exactly -1.
    pub fn is_precursor_im_set(&self) -> bool {
        self.precursor_im != -1.0
    }

    /// The native ID, which is the transition name.
    pub fn native_id(&self) -> &str {
        &self.transition_name
    }

    /// The compound reference (source `getPeptideRef`).
    pub fn peptide_ref(&self) -> &str {
        &self.peptide_ref
    }

    /// The compound reference; the same field as [`Self::peptide_ref`].
    pub fn compound_ref(&self) -> &str {
        &self.peptide_ref
    }

    /// The library intensity.
    pub fn library_intensity(&self) -> f64 {
        self.library_intensity
    }

    /// Set the library intensity.
    pub fn set_library_intensity(&mut self, intensity: f64) {
        self.library_intensity = intensity;
    }

    /// The product m/z.
    pub fn product_mz(&self) -> f64 {
        self.product_mz
    }

    /// The precursor m/z.
    pub fn precursor_mz(&self) -> f64 {
        self.precursor_mz
    }

    /// The precursor ion mobility (-1 when unset).
    pub fn precursor_im(&self) -> f64 {
        self.precursor_im
    }

    /// Set whether the transition is used for detection.
    pub fn set_detecting_transition(&mut self, value: bool) {
        self.flags.set_detecting(value);
    }

    /// Whether the transition is used for detection.
    pub fn is_detecting_transition(&self) -> bool {
        self.flags.detecting()
    }

    /// Set whether the transition is used for quantification.
    pub fn set_quantifying_transition(&mut self, value: bool) {
        self.flags.set_quantifying(value);
    }

    /// Whether the transition is used for quantification.
    pub fn is_quantifying_transition(&self) -> bool {
        self.flags.quantifying()
    }

    /// Set whether the transition is used for identification.
    pub fn set_identifying_transition(&mut self, value: bool) {
        self.flags.set_identifying(value);
    }

    /// Whether the transition is used for identification.
    pub fn is_identifying_transition(&self) -> bool {
        self.flags.identifying()
    }
}

/// A light modification: a UniMod id at a residue location.
///
/// Equality and hashing use both fields.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct LightModification {
    /// Residue location; -1 is the N-terminus and the sequence length the
    /// C-terminus by OpenSWATH convention.
    pub location: i32,
    /// UniMod record id.
    pub unimod_id: i32,
}

/// A light compound: a peptide or a metabolite.
///
/// By convention a compound without a (metabolic) compound name is a peptide.
/// Equality and hashing use [`Self::id`] only (see the module notes).
#[derive(Clone, Debug)]
pub struct LightCompound {
    /// Drift time (ion mobility); -1 means unset.
    pub drift_time: f64,
    /// Retention time in seconds; NaN when unset. The chromatogram extractor
    /// refuses a non-negative extraction window around a NaN `rt`.
    pub rt: f64,
    /// Optional retention-time range start. Both `rt_start` and `rt_end` NaN
    /// means "no range encoded"; otherwise `rt_start..rt_end` is the
    /// per-compound extraction window the chromatogram extractor uses when
    /// called with a NaN extraction window. Filled from a heavyweight peptide
    /// that carries two retention-time entries.
    pub rt_start: f64,
    /// Optional retention-time range end; see [`Self::rt_start`].
    pub rt_end: f64,
    /// Charge; 0 when unset.
    pub charge: i32,
    /// Peptide sequence (empty for a metabolite).
    pub sequence: String,
    /// Ids of the proteins this peptide maps to.
    pub protein_refs: Vec<String>,
    /// Peptide group label (MS:1000893): all peptides that are isotopic forms
    /// of the same peptide should carry the same label.
    pub peptide_group_label: String,
    /// Gene name.
    pub gene_name: String,
    /// Unique compound id, the reference transitions use.
    pub id: String,
    /// Sum formula, for metabolites.
    pub sum_formula: String,
    /// Compound name, for metabolites; empty for peptides.
    pub compound_name: String,
    /// Label type (e.g. `"heavy"` or `"light"`), for round-trip TSV/PQP I/O.
    pub label_type: String,
    /// SMILES representation (metabolomics), for round-trip TSV/PQP I/O.
    pub smiles: String,
    /// Adducts (metabolomics), for round-trip TSV/PQP I/O.
    pub adducts: String,
    /// Modifications of the peptide.
    pub modifications: Vec<LightModification>,
}

impl Default for LightCompound {
    fn default() -> Self {
        Self {
            drift_time: -1.0,
            rt: f64::NAN,
            rt_start: f64::NAN,
            rt_end: f64::NAN,
            charge: 0,
            sequence: String::new(),
            protein_refs: Vec::new(),
            peptide_group_label: String::new(),
            gene_name: String::new(),
            id: String::new(),
            sum_formula: String::new(),
            compound_name: String::new(),
            label_type: String::new(),
            smiles: String::new(),
            adducts: String::new(),
            modifications: Vec::new(),
        }
    }
}

impl PartialEq for LightCompound {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for LightCompound {}

impl Hash for LightCompound {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl LightCompound {
    /// Whether this is a peptide: by convention, a compound without a compound name.
    pub fn is_peptide(&self) -> bool {
        self.compound_name.is_empty()
    }

    /// Set the charge state.
    pub fn set_charge_state(&mut self, charge: i32) {
        self.charge = charge;
    }

    /// The charge state (0 when unset).
    pub fn charge_state(&self) -> i32 {
        self.charge
    }

    /// Set the drift time.
    pub fn set_drift_time(&mut self, drift_time: f64) {
        self.drift_time = drift_time;
    }

    /// The drift time (-1 when unset).
    pub fn drift_time(&self) -> f64 {
        self.drift_time
    }
}

/// A light protein.
///
/// Equality and hashing use [`Self::id`] only (see the module notes).
#[derive(Clone, Debug, Default)]
pub struct LightProtein {
    /// Protein id, the reference compounds use.
    pub id: String,
    /// Protein sequence.
    pub sequence: String,
    /// UniProt identifier, for round-trip TSV/PQP I/O.
    pub uniprot_id: String,
}

impl PartialEq for LightProtein {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for LightProtein {}

impl Hash for LightProtein {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

/// A light targeted experiment: the assay library OpenSWATH extracts and scores.
///
/// # Compound lookup
///
/// The source caches a map from compound id to a raw pointer into
/// `compounds`, built on the first `getCompoundByRef` and never rebuilt: after
/// the compound vector changes (and reallocates), the cached pointers dangle.
/// A reference that no compound carries inserts a null pointer that is then
/// dereferenced. The port keeps no cache inside the value:
/// [`Self::compound_by_ref`] searches, and returns `None` for an unknown
/// reference; [`Self::compound_reference_map`] builds the map once for a
/// caller that looks up many references, borrowing the experiment so it cannot
/// go stale. Both resolve a duplicated id to the **last** compound carrying it,
/// as the source map does (it is filled in index order, later entries
/// overwriting earlier ones).
#[derive(Clone, Debug, Default)]
pub struct LightTargetedExperiment {
    /// All transitions.
    pub transitions: Vec<LightTransition>,
    /// All compounds (peptides and metabolites).
    pub compounds: Vec<LightCompound>,
    /// All proteins.
    pub proteins: Vec<LightProtein>,
}

impl LightTargetedExperiment {
    /// An empty experiment.
    pub fn new() -> Self {
        Self::default()
    }

    /// The transitions.
    pub fn transitions(&self) -> &[LightTransition] {
        &self.transitions
    }

    /// Mutable access to the transitions.
    pub fn transitions_mut(&mut self) -> &mut Vec<LightTransition> {
        &mut self.transitions
    }

    /// The compounds.
    pub fn compounds(&self) -> &[LightCompound] {
        &self.compounds
    }

    /// Mutable access to the compounds.
    pub fn compounds_mut(&mut self) -> &mut Vec<LightCompound> {
        &mut self.compounds
    }

    /// The proteins.
    pub fn proteins(&self) -> &[LightProtein] {
        &self.proteins
    }

    /// Mutable access to the proteins.
    pub fn proteins_mut(&mut self) -> &mut Vec<LightProtein> {
        &mut self.proteins
    }

    /// The compound with id `reference` (the last one when several share it),
    /// or `None`. Linear in the number of compounds; see the type notes.
    pub fn compound_by_ref(&self, reference: &str) -> Option<&LightCompound> {
        self.compounds.iter().rev().find(|c| c.id == reference)
    }

    /// Legacy name of [`Self::compound_by_ref`], source `getPeptideByRef`.
    pub fn peptide_by_ref(&self, reference: &str) -> Option<&LightCompound> {
        self.compound_by_ref(reference)
    }

    /// A map from compound id to compound, the last compound winning for a
    /// duplicated id, for callers that look up many references.
    pub fn compound_reference_map(&self) -> HashMap<&str, &LightCompound> {
        let mut map = HashMap::with_capacity(self.compounds.len());
        for compound in &self.compounds {
            map.insert(compound.id.as_str(), compound);
        }
        map
    }
}

/// A [`LightTransition`] is a kernel [`Transition`](crate::kernel::mrm::Transition),
/// so an OpenSWATH transition group is
/// `MRMTransitionGroup<MSChromatogram, LightTransition>`, as the source's
/// `MRMFeatureFinderScoring::TransitionGroupType` is.
///
/// Native: the source has no such trait; `MRMTransitionGroup` is a template
/// that calls these members by name.
impl crate::kernel::mrm::Transition for LightTransition {
    fn native_id(&self) -> &str {
        &self.transition_name
    }
    fn precursor_mz(&self) -> f64 {
        self.precursor_mz
    }
    fn product_mz(&self) -> f64 {
        self.product_mz
    }
    fn library_intensity(&self) -> f64 {
        self.library_intensity
    }
    fn is_detecting(&self) -> bool {
        self.flags.detecting()
    }
    fn is_identifying(&self) -> bool {
        self.flags.identifying()
    }
    fn is_quantifying(&self) -> bool {
        self.flags.quantifying()
    }
}
