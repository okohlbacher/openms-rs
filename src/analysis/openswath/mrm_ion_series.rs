// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Theoretical fragment ion series for assay generation and decoys.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMIonSeries.h` (with its `.cpp`), the fragment
//! ion bookkeeping shared by `MRMAssay` and `MRMDecoy`. [`MRMIonSeries`]
//! generates theoretical fragment ion series for an [`AASequence`], annotates
//! an observed product m/z against such a series, and annotates a
//! [`ReactionMonitoringTransition`] with the matching PSI-MS CV terms.
//! Neutral losses follow a model similar to the one in SpectraST.
//!
//! # Annotation format
//!
//! The internal annotation format is compatible with SpectraST and derived from
//! Roepstorff and Fohlman (1984, PMID 6525415). An annotation starts with the
//! ion type (`a`, `b`, `c`, `x`, `y`, `z`) followed by the ordinal (the number
//! of residues). A caret `^` followed by a positive integer gives the fragment
//! charge; without a caret the charge is 1. A neutral loss is a `-` followed by
//! either the integer mass (`17` for ammonia) or a composition in
//! [`EmpiricalFormula`] notation (`N1H3` for ammonia).
//!
//! Valid examples: `y3`, `y3^1`, `y3^1-18`, `y3^1-H2O`, `y3-H2O`. As the
//! source's `annotationToCV_` parses them, the charge is split off at the caret
//! first, so the loss belongs before it (`y3-H2O^1`) in anything this module
//! parses; ion-series keys are written in that order.
//!
//! Limitations, as in the source: SpectraST multi-assignments, immonium ions
//! and precursors are not supported.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_ASSAY_SUPPORT.md).
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`MRMIonSeries`]: crate::analysis::openswath::mrm_ion_series::MRMIonSeries
//! [`AASequence`]: crate::chemistry::AASequence
//! [`EmpiricalFormula`]: crate::chemistry::EmpiricalFormula
//! [`ReactionMonitoringTransition`]: crate::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition

use crate::analysis::openswath::reaction_monitoring_transition::ReactionMonitoringTransition;
use crate::analysis::openswath::targeted_experiment_helper::{
    Interpretation, Peptide, UnresolvedModificationPolicy, get_aa_sequence_with,
};
use crate::chemistry::theoretical::residue_loss_formulas;
use crate::chemistry::{AASequence, EmpiricalFormula, ModificationsDB, PeptideFragmentType};
use crate::concept::constants::PROTON_MASS_U;
use crate::concept::math_functions::round_decimal;
use crate::metadata::{CVTerm, MetaValue, MetaValueData};
use crate::param::value::format_float;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

/// An MRM ion series, mapping an annotation (`"y5^2"`, `"b3-H2O1^1"`) to its
/// fragment m/z. The source `IonSeries` is a `std::map`; a [`BTreeMap`] keeps
/// the same key order, which [`MRMIonSeries::annotate_ion`] iterates in.
pub type IonSeries = BTreeMap<String, f64>;

/// The annotation of an ion that matched nothing, and of a transition that
/// could not be annotated.
pub const UNANNOTATED: &str = "unannotated";

/// PSI-MS accession of the `fragment neutral loss` CV term.
pub const FRAGMENT_NEUTRAL_LOSS: &str = "MS:1001524";

/// PSI-MS accession of the `product ion m/z delta` CV term.
pub const PRODUCT_ION_MZ_DELTA: &str = "MS:1000904";

/// Native ceiling on the number of entries [`MRMIonSeries::get_ion_series`]
/// may generate: 2^22. The bound is checked before any entry is computed.
pub const DEFAULT_MAX_ION_SERIES_ENTRIES: usize = 1 << 22;

/// The fragment ion settings the source passes as separate arguments to
/// `getIonSeries`, `annotateTransition`, `MRMAssay` and `MRMDecoy`.
#[derive(Clone, Debug, PartialEq)]
pub struct IonSeriesOptions {
    /// Fragment ion types to generate or accept: any of `a`, `b`, `c`, `x`,
    /// `y`, `z`.
    pub fragment_types: Vec<String>,
    /// Fragment charges to generate or accept. Charges above the precursor
    /// charge are skipped; a charge of zero is rejected.
    pub fragment_charges: Vec<usize>,
    /// Whether specific neutral losses (every declared residue loss other than
    /// the four unspecific ones) are generated.
    pub enable_specific_losses: bool,
    /// Whether the unspecific neutral losses `H2O1`, `H3N1`, `C1H2N2` and
    /// `C1H2N1O1` are generated.
    pub enable_unspecific_losses: bool,
    /// Decimal power that product (and precursor) m/z values are rounded to;
    /// the source default is `-4`.
    pub round_dec_pow: i32,
    /// What converting a TraML peptide does with a modification that cannot be
    /// resolved (see [`get_aa_sequence_with`]). The source skips it with a
    /// warning ([`UnresolvedModificationPolicy::SourceSkip`]); the native
    /// default refuses, because dropping a modification changes the analyte.
    pub unresolved_modifications: UnresolvedModificationPolicy,
    /// Native ceiling on generated ion-series entries, see
    /// [`DEFAULT_MAX_ION_SERIES_ENTRIES`].
    pub max_ion_series_entries: usize,
    /// How formula-style neutral losses in annotations become masses; see
    /// [`FormulaLossParsing`].
    pub formula_loss: FormulaLossParsing,
}

/// How a formula-style neutral loss in an annotation (`y5-H2O1^2`) becomes the
/// mass of its [`FRAGMENT_NEUTRAL_LOSS`] term.
///
/// The source parses the formula into a function-local `static const`
/// (`MRMIonSeries.cpp`, `annotationToCVTermList_`), which C++ initialises only
/// once per process: the first formula loss a process parses fixes the mass
/// reported for every later formula loss, whatever its text. The executed
/// Release build shows it: with `H3N1` parsed first, `y7-H2O1^2` is reported
/// with the ammonia mass -17.0265490957.
#[derive(Clone, Debug, Default)]
pub enum FormulaLossParsing {
    /// Parse every loss formula (native default).
    #[default]
    Parse,
    /// Reproduce the source's `static`: the first formula loss parsed through
    /// this cache fixes the mass of every later one, and later texts are not
    /// parsed (so an invalid one is not an error). A failed first parse leaves
    /// the cache empty, as a throwing `static` initialiser does. Clones share
    /// the cache; use one cache for the scope a source process would span.
    SourceStatic(Arc<Mutex<Option<f64>>>),
}

impl FormulaLossParsing {
    /// A fresh, empty [`FormulaLossParsing::SourceStatic`] cache.
    pub fn source_static() -> Self {
        Self::SourceStatic(Arc::new(Mutex::new(None)))
    }

    /// The loss mass (positive) of the formula `text` under this mode.
    fn mass(&self, text: &str) -> Result<f64> {
        match self {
            Self::Parse => Ok(EmpiricalFormula::parse(text)?.mono_mass()),
            Self::SourceStatic(cache) => {
                let mut guard = cache.lock().unwrap_or_else(|e| e.into_inner());
                match *guard {
                    Some(mass) => Ok(mass),
                    None => {
                        let mass = EmpiricalFormula::parse(text)?.mono_mass();
                        *guard = Some(mass);
                        Ok(mass)
                    }
                }
            }
        }
    }
}

/// Two parse modes are equal when both parse, or both share one cache.
impl PartialEq for FormulaLossParsing {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Parse, Self::Parse) => true,
            (Self::SourceStatic(a), Self::SourceStatic(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }
}

impl IonSeriesOptions {
    /// Options for the given fragment types and charges, with both loss kinds
    /// disabled, the source rounding `-4`, the native
    /// [`UnresolvedModificationPolicy::Reject`] and the default entry ceiling.
    pub fn new<S: Into<String>>(
        fragment_types: impl IntoIterator<Item = S>,
        fragment_charges: impl IntoIterator<Item = usize>,
    ) -> Self {
        Self {
            fragment_types: fragment_types.into_iter().map(Into::into).collect(),
            fragment_charges: fragment_charges.into_iter().collect(),
            enable_specific_losses: false,
            enable_unspecific_losses: false,
            round_dec_pow: -4,
            unresolved_modifications: UnresolvedModificationPolicy::Reject,
            max_ion_series_entries: DEFAULT_MAX_ION_SERIES_ENTRIES,
            formula_loss: FormulaLossParsing::Parse,
        }
    }

    /// The same options with both neutral-loss switches set.
    pub fn with_losses(mut self, specific: bool, unspecific: bool) -> Self {
        self.enable_specific_losses = specific;
        self.enable_unspecific_losses = unspecific;
        self
    }
}

/// Generates theoretical fragment ion series and annotates transitions.
///
/// The source class has no state; its members become associated functions
/// here. See the [module documentation](self) for the annotation format.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MRMIonSeries;

impl MRMIonSeries {
    /// Selects the ion `ionid` from `ionseries`.
    ///
    /// Returns the annotation and its product m/z, or
    /// `("unannotated", -1.0)` when `ionid` is not in the series, as the
    /// source.
    pub fn get_ion(ionseries: &IonSeries, ionid: &str) -> (String, f64) {
        match ionseries.get(ionid) {
            Some(&mz) => (ionid.to_owned(), mz),
            None => (UNANNOTATED.to_owned(), -1.0),
        }
    }

    /// Selects the ion of `ionseries` closest to `product_mz` within
    /// `mz_threshold`.
    ///
    /// Entries are visited in key order. An entry within the threshold
    /// replaces the current best when it is closer by more than 1% of
    /// `mz_threshold`; within that epsilon of the best delta the lower charge
    /// wins, then the lower ordinal (so `b8^2` beats `b18^4`, and `y5` beats
    /// `y15`). A charge that does not parse counts as 1 and an ordinal that
    /// does not parse as `i32::MAX`, as the source helpers.
    ///
    /// Returns the annotation and its product m/z, or `("unannotated", -1.0)`
    /// when nothing is within the threshold. A NaN `product_mz` or threshold
    /// matches nothing.
    pub fn annotate_ion(
        ionseries: &IonSeries,
        product_mz: f64,
        mz_threshold: f64,
    ) -> (String, f64) {
        let mut ion = (UNANNOTATED.to_owned(), -1.0);
        let mut closest_delta = f64::MAX;
        let mut best_charge = i32::MAX;
        let mut best_ordinal = i32::MAX;
        let delta_epsilon = mz_threshold * 0.01;
        for (annotation, &mz) in ionseries {
            let delta = (mz - product_mz).abs();
            if delta > mz_threshold || delta.is_nan() {
                continue;
            }
            let charge = charge_of_annotation(annotation);
            let ordinal = ordinal_of_annotation(annotation);
            let is_better = if delta < closest_delta - delta_epsilon {
                true
            } else if delta <= closest_delta + delta_epsilon {
                charge < best_charge || (charge == best_charge && ordinal < best_ordinal)
            } else {
                false
            };
            if is_better {
                closest_delta = delta;
                best_charge = charge;
                best_ordinal = ordinal;
                ion = (annotation.clone(), mz);
            }
        }
        ion
    }

    /// Annotates `tr` with the CV terms of `annotation`.
    ///
    /// Stores `annotation` as the transition's `annotation` meta value, then
    /// parses it: the product charge is set from the caret (1 without one) and
    /// the product's interpretations are replaced by one with rank 1, the
    /// ordinal, the ion type and, for a neutral loss, a
    /// [`FRAGMENT_NEUTRAL_LOSS`] term holding the negative loss mass (the
    /// integer for a SpectraST loss like `-18`, the monoisotopic formula mass
    /// otherwise). `Precursor_i0` and `MS2_Precursor_i0` yield an empty
    /// interpretation. An unknown ion letter yields
    /// [`PeptideFragmentType::NonIdentified`].
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the charge or ordinal is not an integer
    /// (the source `Exception::ConversionError`), when the ordinal does not fit
    /// the `u8` interpretation ordinal (the source truncates it silently), or
    /// when a loss is neither an integer nor a formula. On error `tr` is
    /// unchanged.
    ///
    /// The source parses a formula loss into a function-local `static`, so a
    /// process caches the *first* formula loss it ever sees and reports that
    /// mass for every later one; this function parses each loss, and
    /// [`annotate_transition_cv_with`](Self::annotate_transition_cv_with)
    /// selects the source behaviour.
    pub fn annotate_transition_cv(
        tr: &mut ReactionMonitoringTransition,
        annotation: &str,
    ) -> Result<()> {
        Self::annotate_transition_cv_with(tr, annotation, &FormulaLossParsing::Parse)
    }

    /// [`annotate_transition_cv`](Self::annotate_transition_cv) with an
    /// explicit [`FormulaLossParsing`] mode.
    ///
    /// # Errors
    ///
    /// As [`annotate_transition_cv`](Self::annotate_transition_cv); under
    /// [`FormulaLossParsing::SourceStatic`] with a filled cache an invalid
    /// loss formula is not an error. On error `tr` is unchanged.
    pub fn annotate_transition_cv_with(
        tr: &mut ReactionMonitoringTransition,
        annotation: &str,
        formula_loss: &FormulaLossParsing,
    ) -> Result<()> {
        let mut next = tr.clone();
        set_annotation(&mut next, annotation);
        annotation_to_cv(&mut next, formula_loss)?;
        *tr = next;
        Ok(())
    }

    /// Annotates the transition `tr` of `peptide` and snaps its m/z values to
    /// theory.
    ///
    /// Without reannotation the transition's first interpretation is used: its
    /// ion type and ordinal select the fragment, whose theoretical m/z is
    /// computed at the product charge (1 when unset); a
    /// [`FRAGMENT_NEUTRAL_LOSS`] term, when either loss switch is on, selects
    /// one of the SpectraST losses -18 (`H2O1`), -17 (`H3N1`), -36 (`H4O2`),
    /// -34 (`H6N2`), -35 (`H5N1O1`), -64 (`C1H4O1S1`), -80 (`H1O3P1`), -98
    /// (`H3O4P1`), -45 (`C1H3N1O1`), -44 (`C1O2`), -46 (`C1H2O2`) or any other
    /// negative mass. The ion type must be in `options.fragment_types` and the
    /// charge in `options.fragment_charges`.
    ///
    /// With `enable_reannotation` the product m/z is instead annotated against
    /// the peptide's full theoretical ion series ([`Self::annotate_ion`]), and
    /// both m/z values are set to their rounded theoretical values.
    ///
    /// The transition is kept annotated only when the product m/z is within
    /// `product_mz_threshold` of theory and the precursor m/z within
    /// `precursor_mz_threshold`; it then gains a [`PRODUCT_ION_MZ_DELTA`] term
    /// and both m/z values are rounded to `options.round_dec_pow`. Otherwise its
    /// `annotation` meta value becomes `"unannotated"` and its interpretation
    /// [`PeptideFragmentType::NonIdentified`].
    ///
    /// # Source behaviour kept
    ///
    /// - The source computes `a`, `b` and `c` ions from the peptide *suffix*
    ///   of the ordinal's length, not the prefix; this is reproduced.
    /// - For a loss outside the eleven named ones the source *adds*
    ///   `|loss| / charge` to the product m/z (`pos -= nl / charge` with a
    ///   negative `nl`), while the named losses subtract; reproduced.
    /// - Reannotation generates the ion series with the default rounding `-4`,
    ///   not `options.round_dec_pow`; reproduced.
    /// - The source's precondition that the peptide has a charge is a debug
    ///   assertion; as its release build, charge 1 is assumed when unset.
    ///
    /// # Errors
    ///
    /// Whatever converting the peptide ([`get_aa_sequence_with`] under
    /// `options.unresolved_modifications`), computing the fragment
    /// ([`AASequence::suffix`], an ordinal past the sequence end) or
    /// [`Self::get_ion_series`] and [`Self::annotate_transition_cv`] return,
    /// and [`Error::InvalidValue`] for a product charge of zero. On error `tr`
    /// is unchanged.
    pub fn annotate_transition(
        tr: &mut ReactionMonitoringTransition,
        peptide: &Peptide,
        precursor_mz_threshold: f64,
        product_mz_threshold: f64,
        enable_reannotation: bool,
        options: &IonSeriesOptions,
    ) -> Result<()> {
        let mut next = tr.clone();
        annotate_transition_into(
            &mut next,
            peptide,
            precursor_mz_threshold,
            product_mz_threshold,
            enable_reannotation,
            options,
        )?;
        *tr = next;
        Ok(())
    }

    /// Computes the theoretical fragment ion series of `sequence`.
    ///
    /// For every fragment type and every charge not above `precursor_charge`,
    /// the ions of ordinals `1..len` are keyed `"<type><ordinal>^<charge>"` with
    /// m/z `(internal residue masses + terminal modification +
    /// charge * proton + ion-type adjustment) / charge`, rounded to
    /// `options.round_dec_pow`.
    /// Prefix ions (`a`, `b`, `c`) carry the N-terminal modification, suffix
    /// ions (`x`, `y`, `z`) the C-terminal one. For every declared residue loss
    /// in the fragment (see the loss table of the support document) an entry
    /// `"<type><ordinal>-<formula>^<charge>"` at `m/z - loss / charge` is added
    /// when the loss is unspecific (`H2O1`, `H3N1`, `C1H2N2`, `C1H2N1O1`) and
    /// `enable_unspecific_losses` is set, or specific and
    /// `enable_specific_losses` is set. An empty sequence yields an empty
    /// series.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] for a fragment type other than `a`, `b`, `c`,
    ///   `x`, `y` or `z` (the source `Exception::IllegalArgument`). As in the
    ///   source it is raised only once that type is reached with a charge not
    ///   above the precursor charge.
    /// - [`Error::InvalidValue`] for a fragment charge of zero, which the
    ///   source divides by.
    /// - [`Error::InvalidValue`] when the series could exceed
    ///   `options.max_ion_series_entries` (native, checked first).
    /// - [`Error::Unsupported`] for a residue without a monoisotopic mass.
    pub fn get_ion_series(
        sequence: &AASequence,
        precursor_charge: usize,
        options: &IonSeriesOptions,
    ) -> Result<IonSeries> {
        let mut ionseries = IonSeries::new();
        let seq_size = sequence.len();
        if seq_size == 0 {
            return Ok(ionseries);
        }
        preflight_ion_series(sequence, precursor_charge, options)?;

        let unspecific = unspecific_losses()?;
        let mut cumulative = Vec::with_capacity(seq_size + 1);
        cumulative.push(0.0_f64);
        let mut losses = Vec::with_capacity(seq_size);
        for index in 0..seq_size {
            let internal = sequence.internal_mass(index)?.ok_or_else(|| {
                Error::Unsupported("residue without a known monoisotopic mass".into())
            })?;
            cumulative.push(cumulative[index] + internal);
            losses.push(residue_loss_formulas(sequence, index)?);
        }
        let n_term = terminal_delta(sequence.n_terminal_modification())?;
        let c_term = terminal_delta(sequence.c_terminal_modification())?;
        let total = cumulative[seq_size];

        for ft in &options.fragment_types {
            for &charge in &options.fragment_charges {
                if charge > precursor_charge {
                    continue;
                }
                let (is_prefix, ion_adj) = match ion_adjustment(ft)? {
                    Some(adjustment) => adjustment,
                    None => {
                        return Err(Error::InvalidValue(format!(
                            "{ft} ion series for peptide sequence \"{sequence}\" with precursor charge +{precursor_charge} could not be generated."
                        )));
                    }
                };
                if charge == 0 {
                    return Err(Error::InvalidValue(
                        "a fragment charge of zero has no m/z".into(),
                    ));
                }
                let charge_f = charge as f64;
                let proton = PROTON_MASS_U * charge_f;
                for i in 1..seq_size {
                    let (internal, range) = if is_prefix {
                        (cumulative[i] + n_term, 0..i)
                    } else {
                        (
                            total - cumulative[seq_size - i] + c_term,
                            seq_size - i..seq_size,
                        )
                    };
                    let pos = (internal + proton + ion_adj) / charge_f;
                    ionseries.insert(
                        format!("{ft}{i}^{charge}"),
                        round_decimal(pos, options.round_dec_pow)?,
                    );
                    for residue_losses in &losses[range] {
                        for loss in residue_losses {
                            let is_unspecific = unspecific.contains(loss);
                            if (options.enable_specific_losses && !is_unspecific)
                                || (options.enable_unspecific_losses && is_unspecific)
                            {
                                ionseries.insert(
                                    format!("{ft}{i}-{loss}^{charge}"),
                                    round_decimal(
                                        pos - loss.mono_mass() / charge_f,
                                        options.round_dec_pow,
                                    )?,
                                );
                            }
                        }
                    }
                }
            }
        }
        Ok(ionseries)
    }
}

/// The four unspecific losses of the source (`H2O1`, `H3N1`, `C1H2N2`,
/// `C1H2N1O1`).
fn unspecific_losses() -> Result<[EmpiricalFormula; 4]> {
    Ok([
        EmpiricalFormula::parse("H2O1")?,
        EmpiricalFormula::parse("H3N1")?,
        EmpiricalFormula::parse("C1H2N2")?,
        EmpiricalFormula::parse("C1H2N1O1")?,
    ])
}

/// Whether `ft` is a prefix type and its `Residue::getInternalTo*Ion` mass;
/// `None` for a type the source rejects.
fn ion_adjustment(ft: &str) -> Result<Option<(bool, f64)>> {
    // Residue.h: NTerm = H, CTerm = OH; a = NTerm - CHO, b = NTerm - H,
    // c = NTerm + NH2, x = CTerm + CO - H, y = CTerm + H, z = CTerm - NH2.
    let (prefix, formula) = match ft {
        "a" => (true, "C-1O-1"),
        "b" => (true, ""),
        "c" => (true, "H3N1"),
        "x" => (false, "C1O2"),
        "y" => (false, "H2O1"),
        "z" => (false, "H-1N-1O1"),
        _ => return Ok(None),
    };
    Ok(Some((
        prefix,
        EmpiricalFormula::parse(formula)?.mono_mass(),
    )))
}

fn terminal_delta(modification: Option<&crate::chemistry::SequenceModification>) -> Result<f64> {
    modification.map_or(Ok(0.0), |m| m.diff_mono_mass())
}

// Upper bound on the entries of one series: every (type, charge) pair yields
// len - 1 ions, each with at most one loss entry per declared residue loss.
fn preflight_ion_series(
    sequence: &AASequence,
    precursor_charge: usize,
    options: &IonSeriesOptions,
) -> Result<()> {
    let overflow = || {
        Error::InvalidValue(format!(
            "ion series would exceed {} entries",
            options.max_ion_series_entries
        ))
    };
    let mut declared_losses = 0usize;
    let losses_enabled = options.enable_specific_losses || options.enable_unspecific_losses;
    for index in (0..sequence.len()).filter(|_| losses_enabled) {
        declared_losses = declared_losses
            .checked_add(residue_loss_formulas(sequence, index)?.len())
            .ok_or_else(overflow)?;
    }
    let charges = options
        .fragment_charges
        .iter()
        .filter(|&&charge| charge <= precursor_charge)
        .count();
    let bound = options
        .fragment_types
        .len()
        .checked_mul(charges)
        .and_then(|n| n.checked_mul(sequence.len() - 1))
        .and_then(|n| n.checked_mul(declared_losses.checked_add(1)?))
        .ok_or_else(overflow)?;
    if bound > options.max_ion_series_entries {
        return Err(overflow());
    }
    Ok(())
}

/// Source `StringUtils::toInt32`:
/// [`crate::concept::string_conversion::to_int32`].
fn to_int32(text: &str) -> Result<i32> {
    crate::concept::string_conversion::to_int32(text).map_err(|_| conversion(text))
}

fn conversion(text: &str) -> Error {
    Error::InvalidValue(format!(
        "Could not convert string '{text}' to an integer value"
    ))
}

/// Charge of an annotation (`b8^2` -> 2, `y5` -> 1); 1 when unparsable.
fn charge_of_annotation(annotation: &str) -> i32 {
    match annotation.find('^') {
        Some(pos) if pos + 1 < annotation.len() => to_int32(&annotation[pos + 1..]).unwrap_or(1),
        _ => 1,
    }
}

/// Ordinal of an annotation (`b8^2` -> 8, `y15` -> 15); `i32::MAX` when
/// absent or unparsable.
fn ordinal_of_annotation(annotation: &str) -> i32 {
    if annotation.is_empty() || annotation == UNANNOTATED {
        return i32::MAX;
    }
    let bytes = annotation.as_bytes();
    let mut end = 1;
    while end < bytes.len() && bytes[end].is_ascii_digit() {
        end += 1;
    }
    if end > 1 {
        to_int32(&annotation[1..end]).unwrap_or(i32::MAX)
    } else {
        i32::MAX
    }
}

fn set_annotation(tr: &mut ReactionMonitoringTransition, annotation: &str) {
    tr.cv_terms
        .metadata
        .insert("annotation".into(), MetaValue::from(annotation));
}

/// The transition's `annotation` meta value as `DataValue::toString` renders
/// it; empty when unset.
fn annotation_of(tr: &ReactionMonitoringTransition) -> String {
    match tr.cv_terms.metadata.get("annotation").map(MetaValue::data) {
        None | Some(MetaValueData::Empty) => String::new(),
        Some(MetaValueData::Float(x)) => format_float(*x, true),
        Some(_) => tr
            .cv_terms
            .metadata
            .get("annotation")
            .map(ToString::to_string)
            .unwrap_or_default(),
    }
}

/// Source `annotationToCV_`: parse the `annotation` meta value into the
/// product charge and a single interpretation.
fn annotation_to_cv(
    tr: &mut ReactionMonitoringTransition,
    formula_loss: &FormulaLossParsing,
) -> Result<()> {
    let text = annotation_of(tr);
    let first = text.split('/').next().unwrap_or_default();
    let mut product = tr.product().clone();
    let annotation = match first.split_once('^') {
        Some((ion, rest)) => {
            // Source splits on every caret and reads the second piece.
            let charge_text = rest.split('^').next().unwrap_or_default();
            product.set_charge_state(to_int32(charge_text)?);
            ion
        }
        None => {
            product.set_charge_state(1);
            first
        }
    };
    let interpretation = annotation_to_interpretation(annotation, formula_loss)?;
    product.reset_interpretations();
    product.add_interpretation(interpretation);
    tr.set_product(product);
    Ok(())
}

/// Source `annotationToCVTermList_`.
fn annotation_to_interpretation(
    annotation: &str,
    formula_loss: &FormulaLossParsing,
) -> Result<Interpretation> {
    let mut interpretation = Interpretation::default();
    let first = annotation.split('/').next().unwrap_or_default();
    if first == "Precursor_i0" || first == "MS2_Precursor_i0" {
        return Ok(interpretation);
    }
    let mut fragment_loss = 0.0;
    let ion = if first.contains('-') {
        let mut pieces = first.split('-');
        let ion = pieces.next().unwrap_or_default();
        let loss = pieces.next().unwrap_or_default();
        // SpectraST integer loss first (boost::lexical_cast<int>: the whole
        // text), else a formula.
        fragment_loss = match loss.parse::<i32>() {
            Ok(nl) => -f64::from(nl),
            Err(_) => -formula_loss.mass(loss)?,
        };
        ion
    } else if first.contains('+') {
        // Neutral gains have no CV term; the source discards them.
        first.split('+').next().unwrap_or_default()
    } else {
        first
    };
    let fragment_type = ion
        .chars()
        .next()
        .ok_or_else(|| Error::InvalidValue("empty fragment annotation".into()))?;
    let ordinal = to_int32(&ion[fragment_type.len_utf8()..])?;
    interpretation.rank = 1;
    interpretation.ordinal = u8::try_from(ordinal).map_err(|_| {
        Error::InvalidValue(format!(
            "fragment ordinal {ordinal} does not fit the interpretation ordinal (0..=255)"
        ))
    })?;
    if fragment_loss < 0.0 {
        let mut term = CVTerm::new(FRAGMENT_NEUTRAL_LOSS, "fragment neutral loss", "MS");
        term.value = MetaValue::try_from(fragment_loss)?;
        interpretation.cv_terms.add(term)?;
    }
    interpretation.iontype = match fragment_type {
        'x' => PeptideFragmentType::XIon,
        'y' => PeptideFragmentType::YIon,
        'z' => PeptideFragmentType::ZIon,
        'a' => PeptideFragmentType::AIon,
        'b' => PeptideFragmentType::BIon,
        'c' => PeptideFragmentType::CIon,
        _ => PeptideFragmentType::NonIdentified,
    };
    Ok(interpretation)
}

/// An `int` charge passed to a `size_t` parameter, as C++ converts it: a
/// negative charge wraps to a huge value (so no fragment charge exceeds it).
#[allow(clippy::cast_sign_loss)]
pub(crate) fn source_size(charge: i32) -> usize {
    charge as isize as usize
}

/// Source `AASequence::getMZ(charge, type)`: the typed monoisotopic ion mass
/// including `charge` protons, divided by the charge.
pub(crate) fn typed_mz(
    sequence: &AASequence,
    fragment: PeptideFragmentType,
    charge: i32,
) -> Result<f64> {
    if charge == 0 {
        return Err(Error::InvalidValue(
            "charge must not be zero for an m/z".into(),
        ));
    }
    Ok(sequence.mono_mass_for(fragment, charge)? / f64::from(charge))
}

/// The value of a CV term as the source reads it back:
/// `StringUtils::toDouble(value.toString())`.
fn cv_value_as_double(term: &CVTerm) -> Result<f64> {
    let text = match term.value.data() {
        MetaValueData::Float(x) => format_float(*x, true),
        _ => term.value.to_string(),
    };
    text.trim_matches([' ', '\t', '\n', '\r'])
        .parse::<f64>()
        .map_err(|_| Error::InvalidValue(format!("could not convert '{text}' to a double")))
}

// The eleven SpectraST losses the source names, by integer mass.
const NAMED_LOSSES: [(f64, &str); 11] = [
    (-18.0, "H2O1"),
    (-17.0, "H3N1"),
    (-36.0, "H4O2"),
    (-34.0, "H6N2"),
    (-35.0, "H5N1O1"),
    (-64.0, "C1H4O1S1"),
    (-80.0, "H1O3P1"),
    (-98.0, "H3O4P1"),
    (-45.0, "C1H3N1O1"),
    (-44.0, "C1O2"),
    (-46.0, "C1H2O2"),
];

fn annotate_transition_into(
    tr: &mut ReactionMonitoringTransition,
    peptide: &Peptide,
    precursor_mz_threshold: f64,
    product_mz_threshold: f64,
    enable_reannotation: bool,
    options: &IonSeriesOptions,
) -> Result<()> {
    let sequence = get_aa_sequence_with(
        peptide,
        ModificationsDB::global(),
        options.unresolved_modifications,
    )?
    .sequence;
    let precursor_charge = peptide.charge_state().unwrap_or(1);
    let fragment_charge = tr.product_charge_state().unwrap_or(1);
    let mut prec_pos = typed_mz(&sequence, PeptideFragmentType::Full, precursor_charge)?;
    let mut unannotated = false;
    let mut pos = -1.0;
    let mut ionstring = String::new();
    let mut interpretation = Interpretation::default();

    if let Some(first) = tr.product().interpretation_list().first() {
        interpretation = first.clone();
        if interpretation.ordinal > 0 {
            let ordinal = usize::from(interpretation.ordinal);
            let letter = match interpretation.iontype {
                PeptideFragmentType::XIon => Some("x"),
                PeptideFragmentType::YIon => Some("y"),
                PeptideFragmentType::ZIon => Some("z"),
                PeptideFragmentType::AIon => Some("a"),
                PeptideFragmentType::BIon => Some("b"),
                PeptideFragmentType::CIon => Some("c"),
                _ => None,
            };
            match letter {
                Some(letter) => {
                    // Source: getSuffix for every ion type, a/b/c included.
                    let ion = sequence.suffix(ordinal)?;
                    ionstring.push_str(letter);
                    pos = typed_mz(&ion, interpretation.iontype, fragment_charge)?;
                }
                None => unannotated = true,
            }
        } else {
            unannotated = true;
        }
        if !options.fragment_types.contains(&ionstring) {
            unannotated = true;
        }
        if interpretation.ordinal > 0 {
            ionstring.push_str(&interpretation.ordinal.to_string());
        } else {
            unannotated = true;
        }
        let loss_term = interpretation
            .cv_terms
            .get(FRAGMENT_NEUTRAL_LOSS)
            .and_then(|terms| terms.first());
        if let Some(term) =
            loss_term.filter(|_| options.enable_specific_losses || options.enable_unspecific_losses)
        {
            let nl = cv_value_as_double(term)?;
            let charge = f64::from(fragment_charge);
            if let Some((_, formula)) = NAMED_LOSSES.iter().find(|(mass, _)| *mass == nl) {
                ionstring.push('-');
                ionstring.push_str(formula);
                pos -= EmpiricalFormula::parse(formula)?.mono_mass() / charge;
            } else if nl < 0.0 {
                // Source: toStr(roundDecimal(nl)) supplies the minus sign, and
                // `pos -= nl / charge` adds the loss mass back.
                ionstring.push_str(&format_float(
                    round_decimal(nl, options.round_dec_pow)?,
                    true,
                ));
                pos -= nl / charge;
            } else {
                unannotated = true;
            }
        }
        let charge_accepted = usize::try_from(fragment_charge)
            .ok()
            .filter(|&c| c >= 1)
            .is_some_and(|c| options.fragment_charges.contains(&c));
        if charge_accepted {
            ionstring.push('^');
            ionstring.push_str(&fragment_charge.to_string());
            set_annotation(tr, &ionstring);
        } else {
            unannotated = true;
        }
    } else {
        unannotated = true;
    }

    if enable_reannotation {
        let mut series_options = options.clone();
        // Source: getIonSeries is called without round_decPow here.
        series_options.round_dec_pow = -4;
        let precursor = source_size(precursor_charge);
        let ionseries = MRMIonSeries::get_ion_series(&sequence, precursor, &series_options)?;
        let (annotation, mz) =
            MRMIonSeries::annotate_ion(&ionseries, tr.product_mz(), product_mz_threshold);
        ionstring = annotation;
        set_annotation(tr, &ionstring);
        pos = round_decimal(mz, options.round_dec_pow)?;
        prec_pos = round_decimal(prec_pos, options.round_dec_pow)?;
        tr.set_product_mz(pos);
        tr.set_precursor_mz(prec_pos);
        if ionstring == UNANNOTATED {
            unannotated = true;
        } else {
            annotation_to_cv(tr, &options.formula_loss)?;
            interpretation = tr
                .product()
                .interpretation_list()
                .first()
                .cloned()
                .unwrap_or_default();
            unannotated = false;
        }
    }

    if !unannotated
        && (tr.product_mz() - pos).abs() <= product_mz_threshold
        && (tr.precursor_mz() - prec_pos).abs() <= precursor_mz_threshold
    {
        let mut term = CVTerm::new(PRODUCT_ION_MZ_DELTA, "product ion m/z delta", "MS");
        term.value = MetaValue::try_from(
            round_decimal(tr.product_mz() - pos, options.round_dec_pow)?.abs(),
        )?;
        interpretation.cv_terms.replace(term)?;
        pos = round_decimal(pos, options.round_dec_pow)?;
        prec_pos = round_decimal(prec_pos, options.round_dec_pow)?;
        tr.set_product_mz(pos);
        tr.set_precursor_mz(prec_pos);
    } else {
        unannotated = true;
    }

    if unannotated {
        interpretation.iontype = PeptideFragmentType::NonIdentified;
        set_annotation(tr, UNANNOTATED);
    } else {
        set_annotation(tr, &ionstring);
        annotation_to_cv(tr, &options.formula_loss)?;
    }
    let mut product = tr.product().clone();
    product.reset_interpretations();
    product.add_interpretation(interpretation);
    tr.set_product(product);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn to_int32_follows_the_source_conversion() {
        assert_eq!(to_int32(" +12 ").unwrap(), 12);
        assert_eq!(to_int32("-3").unwrap(), -3);
        assert!(to_int32("").is_err());
        assert!(to_int32("++1").is_err());
        assert!(to_int32("1a").is_err());
    }

    #[test]
    fn annotation_helpers_match_the_source_fallbacks() {
        assert_eq!(charge_of_annotation("b8^2"), 2);
        assert_eq!(charge_of_annotation("y5"), 1);
        assert_eq!(charge_of_annotation("y5^x"), 1);
        assert_eq!(ordinal_of_annotation("y15^3"), 15);
        assert_eq!(ordinal_of_annotation("unannotated"), i32::MAX);
        assert_eq!(ordinal_of_annotation("y"), i32::MAX);
    }
}
