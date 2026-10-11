// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Assay generation from a raw targeted experiment.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMAssay.h` (with its `.cpp`), the algorithm of
//! the `OpenSwathAssayGenerator` TOPP tool (and the compound filters of
//! `TargetedFileConverter`-style metabolomics workflows). [`MRMAssay`]
//! generates assays from a raw, unfiltered [`TargetedExperiment`] as
//! `TargetedFileConverter` produces it.
//!
//! Transitions are selected by a set of rules described in Schubert et al.,
//! 2015 (PMID 25675208): [`MRMAssay::reannotate_transitions`] annotates them
//! against theory, [`MRMAssay::restrict_transitions`] drops unannotated ones and
//! those outside an m/z range or inside the precursor's isolation window, and
//! [`MRMAssay::detecting_transitions`] keeps the most intense ones per
//! precursor. In addition, [`MRMAssay::uis_transitions`] generates unique ion
//! signature (UIS) transitions (Sherman et al., 2009, PMID 19556279) from in
//! silico ion series, as described for IPF in Rosenberger et al., 2017 (PMID
//! 28604659).
//!
//! # Not ported
//!
//! The source's `LightTargetedExperiment` methods (`reannotateTransitionsLight`,
//! `restrictTransitionsLight`, `detectingTransitionsLight`,
//! `uisTransitionsLight` and their helpers) operate on
//! `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`, which is not ported yet.
//! The `ProgressLogger` base is not ported; progress is not reported. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_ASSAY_SUPPORT.md).
//!
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`MRMAssay`]: crate::analysis::openswath::mrm_assay::MRMAssay
//! [`MRMAssay::reannotate_transitions`]: crate::analysis::openswath::mrm_assay::MRMAssay::reannotate_transitions
//! [`MRMAssay::restrict_transitions`]: crate::analysis::openswath::mrm_assay::MRMAssay::restrict_transitions
//! [`MRMAssay::detecting_transitions`]: crate::analysis::openswath::mrm_assay::MRMAssay::detecting_transitions
//! [`MRMAssay::uis_transitions`]: crate::analysis::openswath::mrm_assay::MRMAssay::uis_transitions
//! [`TargetedExperiment`]: crate::analysis::openswath::targeted_experiment::TargetedExperiment

use crate::analysis::openswath::mrm_decoy::{BoostUniformInt, RANDOM_RESIDUES, time_seed};
use crate::analysis::openswath::mrm_ion_series::{
    IonSeriesOptions, MRMIonSeries, UNANNOTATED, source_size, typed_mz,
};
use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{
    Compound, Peptide, Protein, get_aa_sequence_with,
};
use crate::chemistry::{
    AASequence, ModificationsDB, PeptideFragmentType, ResidueModification, SequenceModification,
    TermSpecificity,
};
use crate::concept::math_functions::round_decimal;
use crate::metadata::MetaValue;
use crate::param::value::format_float;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashSet};

/// A fragment sequence map: `(fragment m/z, modified sequence)` pairs (source
/// `FragmentSeqMap`).
pub type FragmentSeqMap = Vec<(f64, String)>;

/// An ion series as a list: `(annotation, fragment m/z)` pairs (source
/// `MRMAssay::IonSeries`).
pub type AssayIonSeries = Vec<(String, f64)>;

/// Native ceiling on the n-choose-k combinations one modification may
/// enumerate while generating peptidoforms: 2^20. The source enumerates them
/// all before comparing against `max_num_alternative_localizations`.
pub const MAX_PEPTIDOFORM_COMBINATIONS: usize = 1 << 20;

/// Native ceiling on the transitions [`MRMAssay::uis_transitions`] may
/// append: 2^26.
pub const MAX_UIS_TRANSITIONS: usize = 1 << 26;

/// The UIS settings of the source `uisTransitions` besides the ion series.
#[derive(Clone, Debug, PartialEq)]
pub struct UisOptions {
    /// Fragment types, charges, losses, rounding and modification policy.
    /// As the source, the in silico ion series are generated with rounding
    /// `-4` and then rounded again to `ion_series.round_dec_pow`.
    pub ion_series: IonSeriesOptions,
    /// Add the precursor (`MS2_Precursor_i0`) as a theoretical transition.
    pub enable_ms2_precursors: bool,
    /// Product m/z threshold in Th within which two ions interfere.
    pub mz_threshold: f64,
    /// SWATH isolation windows `(lower, upper)`, sorted by lower bound; empty
    /// for none.
    pub swathes: Vec<(f64, f64)>,
    /// Maximum number of alternative peptidoforms; peptides with more are
    /// skipped. Source default 20.
    pub max_num_alternative_localizations: usize,
    /// Seed of the random decoy sequences. `None` seeds from the current time
    /// in seconds, as the source's `-1` (its default).
    pub shuffle_seed: Option<u32>,
    /// Skip the generation of decoy UIS transitions.
    pub disable_decoy_transitions: bool,
}

impl UisOptions {
    /// The source defaults for everything but the ion series and the m/z
    /// threshold: no MS2 precursors, no swathes, 20 alternative
    /// localizations, time-based seed, decoys enabled.
    pub fn new(ion_series: IonSeriesOptions, mz_threshold: f64) -> Self {
        Self {
            ion_series,
            enable_ms2_precursors: false,
            mz_threshold,
            swathes: Vec::new(),
            max_num_alternative_localizations: 20,
            shuffle_seed: None,
            disable_decoy_transitions: false,
        }
    }
}

/// Generates and filters assays of a [`TargetedExperiment`].
///
/// The source class has no state besides its `ProgressLogger` base; its
/// members become associated functions here.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MRMAssay;

// Keys of the source's `SequenceMapT` and `IonMapT`: the `int` swath index
// stored as `size_t`, so "no swath" (-1) sorts last.
type SwathKey = usize;
type SequenceMap = BTreeMap<SwathKey, BTreeMap<String, BTreeSet<String>>>;
type IonMap = BTreeMap<SwathKey, BTreeMap<String, FragmentSeqMap>>;
type PeptideMap = BTreeMap<String, AssayIonSeries>;

impl MRMAssay {
    /// The peptidoforms whose fragments lie within `mz_threshold` of
    /// `fragment_ion` (source `getMatchingPeptidoforms_`).
    ///
    /// `ions` must be sorted by m/z ascending; the search is a binary search
    /// for the lower bound followed by a scan. The result is sorted and free
    /// of duplicates.
    pub fn matching_peptidoforms(
        fragment_ion: f64,
        ions: &[(f64, String)],
        mz_threshold: f64,
    ) -> Vec<String> {
        let lower = fragment_ion - mz_threshold;
        let upper = fragment_ion + mz_threshold;
        let start = ions.partition_point(|(mz, _)| *mz < lower);
        let mut isoforms: Vec<String> = ions[start..]
            .iter()
            .take_while(|(mz, _)| *mz <= upper)
            .map(|(_, sequence)| sequence.clone())
            .collect();
        isoforms.sort();
        isoforms.dedup();
        isoforms
    }

    /// The index of the swath containing `precursor_mz`, or `-1` (source
    /// `getSwath_`).
    ///
    /// `swathes` must be sorted by lower bound. Bounds are inclusive; in an
    /// overlap the later (upper) swath wins.
    pub fn swath_index(swathes: &[(f64, f64)], precursor_mz: f64) -> i32 {
        let upper = swathes
            .partition_point(|(lower, _)| precursor_mz.partial_cmp(lower) != Some(Ordering::Less));
        let contains = |i: usize| {
            swathes
                .get(i)
                .is_some_and(|(lo, hi)| precursor_mz >= *lo && precursor_mz <= *hi)
        };
        let mut swath = -1;
        if upper > 0 && contains(upper - 1) {
            swath = i32::try_from(upper - 1).unwrap_or(-1);
        }
        if contains(upper) {
            swath = i32::try_from(upper).unwrap_or(-1);
        }
        swath
    }

    /// Whether `product_mz` falls into the swath of `precursor_mz` (source
    /// `isInSwath_`). A precursor outside every swath counts as inside, so the
    /// transition is removed.
    pub fn is_in_swath(swathes: &[(f64, f64)], precursor_mz: f64, product_mz: f64) -> bool {
        let index = Self::swath_index(swathes, precursor_mz);
        match usize::try_from(index).ok().and_then(|i| swathes.get(i)) {
            None => true,
            Some((lo, hi)) => product_mz >= *lo && product_mz <= *hi,
        }
    }

    /// A random sequence of `sequence_size` residues (source
    /// `getRandomSequence_`), each `ANDCEQGHILMFSTWYV[draw % 17]` with the
    /// source's `[0, 9]` draw, so only the first ten residues occur.
    pub(crate) fn random_sequence(sequence_size: usize, rng: &mut BoostUniformInt) -> String {
        (0..sequence_size)
            .map(|_| char::from(RANDOM_RESIDUES[(rng.draw() % 17) as usize]))
            .collect()
    }

    /// A random sequence of `sequence_size` residues (source
    /// `getRandomSequence_`) from a fresh `boost::mt19937(seed)`, as the class
    /// test drives it: each residue is `ANDCEQGHILMFSTWYV[draw % 17]` with the
    /// source's `[0, 9]` draw, so only the first ten residues occur.
    pub fn random_sequence_with_seed(sequence_size: usize, seed: u32) -> String {
        Self::random_sequence(sequence_size, &mut BoostUniformInt::new(seed))
    }

    /// All `k`-element combinations of `n`, in the source's order (source
    /// `nchoosekcombinations_`: a bitmask of `k` ones stepped with
    /// `prev_permutation`).
    ///
    /// As the source, `k > n.len()` yields the single combination of all of
    /// `n`, and `k == 0` the single empty combination.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when more than [`MAX_PEPTIDOFORM_COMBINATIONS`]
    /// combinations would result (native).
    pub fn n_choose_k_combinations(n: &[usize], k: usize) -> Result<Vec<Vec<usize>>> {
        let k = k.min(n.len());
        if binomial_exceeds(n.len(), k, MAX_PEPTIDOFORM_COMBINATIONS) {
            return Err(Error::InvalidValue(format!(
                "more than {MAX_PEPTIDOFORM_COMBINATIONS} modification combinations"
            )));
        }
        let mut bitmask: Vec<u8> = (0..n.len()).map(|i| u8::from(i < k)).collect();
        let mut combinations = Vec::new();
        loop {
            combinations.push(
                n.iter()
                    .zip(&bitmask)
                    .filter(|(_, bit)| **bit == 1)
                    .map(|(&value, _)| value)
                    .collect(),
            );
            if !prev_permutation(&mut bitmask) {
                break;
            }
        }
        Ok(combinations)
    }

    /// Apply `modification` at every combination of `mods_combs` to every
    /// sequence (source `addModificationsSequences_`).
    ///
    /// Positions are 1-based residues; `0` is the N-terminus and `len + 1`
    /// the C-terminus. A terminal position is applied only when the first
    /// terminal record of `modification` in the registry names the terminal
    /// residue or `X`; otherwise that combination is dropped. A combination
    /// that would modify an already modified residue is dropped.
    ///
    /// # Errors
    ///
    /// Whatever setting the modification on an [`AASequence`] returns (an
    /// unknown name or incompatible residue), and [`Error::InvalidValue`] for
    /// a position past the C-terminus.
    pub fn add_modifications_sequences(
        sequences: &[AASequence],
        mods_combs: &[Vec<usize>],
        modification: &str,
    ) -> Result<Vec<AASequence>> {
        let db = ModificationsDB::global();
        let nterm = db.find(modification, None, Some(TermSpecificity::NTerm));
        let cterm = db.find(modification, None, Some(TermSpecificity::CTerm));
        let origin_matches = |records: &[&ResidueModification], residue: Option<char>| {
            records.first().is_some_and(|m| match m.origin() {
                None | Some('X') => true,
                Some(origin) => residue == Some(origin),
            })
        };
        let mut modified = Vec::new();
        for sequence in sequences {
            for combination in mods_combs {
                let mut multi_mod = false;
                let mut skip_invalid = false;
                let mut temp = sequence.clone();
                for &position in combination {
                    let len = temp.len();
                    if position == 0 {
                        if origin_matches(&nterm, temp.as_str().chars().next()) {
                            temp.set_n_terminal_modification(modification)?;
                        } else {
                            skip_invalid = true;
                        }
                    } else if position == len + 1 {
                        if origin_matches(&cterm, temp.as_str().chars().last()) {
                            temp.set_c_terminal_modification(modification)?;
                        } else {
                            skip_invalid = true;
                        }
                    } else if position > len + 1 {
                        return Err(Error::InvalidValue(format!(
                            "modification position {position} is past the C-terminus"
                        )));
                    } else if temp.residue_modification(position - 1)?.is_none() {
                        temp.set_modification(position - 1, modification)?;
                    } else {
                        multi_mod = true;
                    }
                }
                if !skip_invalid && !multi_mod {
                    modified.push(temp);
                }
            }
        }
        Ok(modified)
    }

    /// All alternative peptidoforms of `sequence` compatible with the
    /// modification registry (source `generateTheoreticalPeptidoforms_`).
    ///
    /// Every modification of `sequence` (by name, with its count) is placed
    /// at every n-choose-k combination of the positions that can carry it:
    /// the N-terminus when it has an N-terminal record, the C-terminus when it
    /// has a C-terminal one, and every residue with a matching
    /// residue-specific record. Modifications are applied in name order.
    ///
    /// # Errors
    ///
    /// As [`n_choose_k_combinations`](Self::n_choose_k_combinations) and
    /// [`add_modifications_sequences`](Self::add_modifications_sequences).
    pub fn theoretical_peptidoforms(sequence: &AASequence) -> Result<Vec<AASequence>> {
        Self::theoretical_peptidoforms_decoy(sequence, sequence)
    }

    /// As [`theoretical_peptidoforms`](Self::theoretical_peptidoforms), but
    /// the permutations inferred from `sequence` are applied to the unmodified
    /// `decoy_sequence`, which may have more modifiable residues (source
    /// `generateTheoreticalPeptidoformsDecoy_`). For example the target
    /// `SAS(Phospho)K` gives `[SAS(Phospho)K, S(Phospho)ASK]`, and the decoy
    /// `SSS(Phospho)K` then gives `[SSS(Phospho)K, S(Phospho)SSK]`.
    ///
    /// # Errors
    ///
    /// As [`theoretical_peptidoforms`](Self::theoretical_peptidoforms).
    pub fn theoretical_peptidoforms_decoy(
        sequence: &AASequence,
        decoy_sequence: &AASequence,
    ) -> Result<Vec<AASequence>> {
        let db = ModificationsDB::global();
        let mut mods: BTreeMap<String, usize> = BTreeMap::new();
        for modification in [
            sequence.n_terminal_modification(),
            sequence.c_terminal_modification(),
        ]
        .into_iter()
        .flatten()
        {
            *mods.entry(modification_name(modification)).or_default() += 1;
        }
        for i in 0..sequence.len() {
            if let Some(modification) = sequence.residue_modification(i)? {
                *mods.entry(modification_name(modification)).or_default() += 1;
            }
        }
        let mut sequences = vec![AASequence::parse(decoy_sequence.as_str())?];
        for (name, &count) in &mods {
            let mut positions = Vec::new();
            if !db.find(name, None, Some(TermSpecificity::NTerm)).is_empty() {
                positions.push(0);
            }
            if !db.find(name, None, Some(TermSpecificity::CTerm)).is_empty() {
                positions.push(sequence.len() + 1);
            }
            for (i, residue) in sequence.as_str().chars().enumerate() {
                if !db
                    .find(name, Some(residue), Some(TermSpecificity::Anywhere))
                    .is_empty()
                {
                    positions.push(i + 1);
                }
            }
            let combinations = Self::n_choose_k_combinations(&positions, count)?;
            sequences = Self::add_modifications_sequences(&sequences, &combinations, name)?;
        }
        Ok(sequences)
    }

    /// Annotate the transitions of `exp` against the theoretical ion series of
    /// their peptides and keep only the annotated ones (source
    /// `reannotateTransitions`).
    ///
    /// Transitions are processed grouped by peptide reference in ascending
    /// order. Each is annotated with [`MRMIonSeries::annotate_ion`] within
    /// `product_mz_threshold`; a precursor m/z further than
    /// `precursor_mz_threshold` from the rounded theoretical value makes it
    /// unannotated. Kept transitions get the theoretical (rounded) precursor
    /// and product m/z, the CV terms of their annotation
    /// ([`MRMIonSeries::annotate_transition_cv`]) and their peptide reference.
    /// Only the transition list of `exp` changes.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when a transition's peptide is missing,
    /// and whatever peptide conversion and the ion series return. On error
    /// `exp` is unchanged.
    pub fn reannotate_transitions(
        exp: &mut TargetedExperiment,
        precursor_mz_threshold: f64,
        product_mz_threshold: f64,
        options: &IonSeriesOptions,
    ) -> Result<()> {
        let db = ModificationsDB::global();
        let mut transitions = Vec::new();
        for (peptide_ref, group) in group_by(exp.transitions(), |t| t.peptide_ref()) {
            let peptide = exp.peptide_by_ref(peptide_ref)?;
            let sequence =
                get_aa_sequence_with(peptide, db, options.unresolved_modifications)?.sequence;
            let charge = peptide.charge_state().unwrap_or(1);
            let series = MRMIonSeries::get_ion_series(&sequence, source_size(charge), options)?;
            let precursor_mz = round_decimal(
                typed_mz(&sequence, PeptideFragmentType::Full, charge)?,
                options.round_dec_pow,
            )?;
            for tr in group {
                let (mut annotation, mz) =
                    MRMIonSeries::annotate_ion(&series, tr.product_mz(), product_mz_threshold);
                if (tr.precursor_mz() - precursor_mz).abs() > precursor_mz_threshold {
                    annotation = UNANNOTATED.into();
                }
                if annotation == UNANNOTATED {
                    continue;
                }
                let mut tr = tr.clone();
                tr.set_precursor_mz(precursor_mz);
                tr.set_product_mz(mz);
                MRMIonSeries::annotate_transition_cv_with(
                    &mut tr,
                    &annotation,
                    &options.formula_loss,
                )?;
                tr.set_peptide_ref(peptide.id.clone());
                transitions.push(tr);
            }
        }
        exp.set_transitions(transitions);
        Ok(())
    }

    /// Drop unannotated transitions, those whose product m/z falls into the
    /// precursor's swath, and those outside `[lower_mz_limit,
    /// upper_mz_limit]` (source `restrictTransitions`).
    ///
    /// A transition is unannotated when its first interpretation has type
    /// [`PeptideFragmentType::NonIdentified`]; one without interpretations is
    /// kept by that rule. With `swathes` empty the swath rule is off; with
    /// swathes, a precursor outside all of them removes the transition.
    ///
    /// # Errors
    ///
    /// [`Error::MissingInformation`] when a transition's peptide is missing,
    /// and whatever converting it returns: the source converts every peptide
    /// (for its debug log) and so fails the same way. On error `exp` is
    /// unchanged.
    pub fn restrict_transitions(
        exp: &mut TargetedExperiment,
        lower_mz_limit: f64,
        upper_mz_limit: f64,
        swathes: &[(f64, f64)],
        options: &IonSeriesOptions,
    ) -> Result<()> {
        let db = ModificationsDB::global();
        let mut transitions = Vec::with_capacity(exp.transitions().len());
        for tr in exp.transitions() {
            let peptide = exp.peptide_by_ref(tr.peptide_ref())?;
            get_aa_sequence_with(peptide, db, options.unresolved_modifications)?;
            if tr
                .product()
                .interpretation_list()
                .first()
                .is_some_and(|i| i.iontype == PeptideFragmentType::NonIdentified)
            {
                continue;
            }
            if !swathes.is_empty() && Self::is_in_swath(swathes, tr.precursor_mz(), tr.product_mz())
            {
                continue;
            }
            if tr.product_mz() < lower_mz_limit || tr.product_mz() > upper_mz_limit {
                continue;
            }
            transitions.push(tr.clone());
        }
        exp.set_transitions(transitions);
        Ok(())
    }

    /// Select the detecting transitions: per precursor with at least
    /// `min_transitions` transitions, the (at most `max_transitions`) most
    /// intense non-decoy ones (source `detectingTransitions`).
    ///
    /// Transitions are grouped by peptide reference in ascending order. The
    /// library intensities of a group are sorted descending and cut to
    /// `max_transitions`; a transition whose intensity is among them is kept
    /// and marked detecting, in input order, until `max_transitions` are kept
    /// (so ties do not exceed the limit). Peptides without kept transitions
    /// and proteins without kept peptides are removed.
    ///
    /// The source's `int` limits are compared as `size_t`, so a negative
    /// `min_transitions` keeps nothing and a negative `max_transitions` no
    /// limit; reproduced.
    pub fn detecting_transitions(
        exp: &mut TargetedExperiment,
        min_transitions: i32,
        max_transitions: i32,
    ) {
        let (min, max) = (source_size(min_transitions), source_size(max_transitions));
        let mut transitions = Vec::new();
        let mut peptide_ids: HashSet<String> = HashSet::new();
        for (_, group) in group_by(exp.transitions(), |t| t.peptide_ref()) {
            if group.len() < min {
                continue;
            }
            let top = top_intensities(&group, max);
            let mut kept = 0usize;
            for tr in group {
                if top.contains(&tr.library_intensity().to_bits())
                    && tr.decoy_transition_type() != DecoyTransitionType::Decoy
                    && kept < max
                {
                    let mut tr = tr.clone();
                    tr.set_detecting_transition(true);
                    kept += 1;
                    peptide_ids.insert(tr.peptide_ref().to_owned());
                    transitions.push(tr);
                }
            }
        }
        let mut protein_ids: HashSet<String> = HashSet::new();
        let peptides: Vec<Peptide> = exp
            .peptides()
            .iter()
            .filter(|p| peptide_ids.contains(&p.id))
            .inspect(|p| protein_ids.extend(p.protein_refs.iter().cloned()))
            .cloned()
            .collect();
        let proteins: Vec<Protein> = exp
            .proteins()
            .iter()
            .filter(|p| protein_ids.contains(&p.id))
            .cloned()
            .collect();
        exp.set_transitions(transitions);
        exp.set_peptides(peptides);
        exp.set_proteins(proteins);
    }

    /// Keep the (at most `max_transitions`) most intense transitions per
    /// compound with at least `min_transitions` transitions, or whose first
    /// transition is a decoy (source `filterMinMaxTransitionsCompound`).
    ///
    /// Grouping, intensity ranking and the `int` limits are as in
    /// [`detecting_transitions`](Self::detecting_transitions), by compound
    /// reference; decoys are not excluded from the ranking. Compounds without
    /// kept transitions are removed.
    pub fn filter_min_max_transitions_compound(
        exp: &mut TargetedExperiment,
        min_transitions: i32,
        max_transitions: i32,
    ) {
        let (min, max) = (source_size(min_transitions), source_size(max_transitions));
        let mut transitions = Vec::new();
        let mut compound_ids: Vec<String> = Vec::new();
        for (_, group) in group_by(exp.transitions(), |t| t.compound_ref()) {
            let first_is_decoy = group
                .first()
                .is_some_and(|t| t.decoy_transition_type() == DecoyTransitionType::Decoy);
            if !(group.len() >= min || first_is_decoy) {
                continue;
            }
            let top = top_intensities(&group, max);
            let mut kept = 0usize;
            for tr in group {
                if top.contains(&tr.library_intensity().to_bits()) && kept < max {
                    let mut tr = tr.clone();
                    tr.set_detecting_transition(true);
                    kept += 1;
                    if !compound_ids.iter().any(|id| id == tr.compound_ref()) {
                        compound_ids.push(tr.compound_ref().to_owned());
                    }
                    transitions.push(tr);
                }
            }
        }
        let compounds: Vec<Compound> = exp
            .compounds()
            .iter()
            .filter(|c| compound_ids.iter().any(|id| *id == c.id))
            .cloned()
            .collect();
        exp.set_transitions(transitions);
        exp.set_compounds(compounds);
    }

    /// Remove decoy compounds (and their transitions) that have no target
    /// (source `filterUnreferencedDecoysCompound`).
    ///
    /// A compound whose id contains `decoy` is a decoy; its target id is the
    /// id with the first `_decoy` removed, so
    /// `84_CompoundName_decoy_[M+H]+_88_22` references
    /// `84_CompoundName_[M+H]+_88_22`. Decoys whose target id is not a
    /// compound id are removed, as are the transitions referencing them.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when an id contains `decoy` but not `_decoy`,
    /// where the source's `std::string::erase(npos, ...)` throws
    /// `std::out_of_range`. On error `exp` is unchanged.
    pub fn filter_unreferenced_decoys_compound(exp: &mut TargetedExperiment) -> Result<()> {
        const DECOY_SUFFIX: &str = "_decoy";
        let mut targets: Vec<String> = Vec::new();
        let mut decoys: Vec<String> = Vec::new();
        let mut references: Vec<(String, String)> = Vec::new();
        for compound in exp.compounds() {
            if compound.id.contains("decoy") {
                let position = compound.id.find(DECOY_SUFFIX).ok_or_else(|| {
                    Error::InvalidValue(format!(
                        "decoy compound id '{}' does not contain '{DECOY_SUFFIX}'",
                        compound.id
                    ))
                })?;
                let mut target = compound.id.clone();
                target.replace_range(position..position + DECOY_SUFFIX.len(), "");
                decoys.push(target.clone());
                references.push((compound.id.clone(), target));
            } else {
                targets.push(compound.id.clone());
            }
        }
        targets.sort();
        decoys.sort();
        let missing = multiset_difference(&decoys, &targets);
        let single_decoy_ids: HashSet<&str> = missing
            .iter()
            .filter_map(|target| {
                references
                    .iter()
                    .find(|(_, t)| t == target)
                    .map(|(decoy, _)| decoy.as_str())
            })
            .collect();
        let compounds: Vec<Compound> = exp
            .compounds()
            .iter()
            .filter(|c| !single_decoy_ids.contains(c.id.as_str()))
            .cloned()
            .collect();
        let transitions: Vec<ReactionMonitoringTransition> = exp
            .transitions()
            .iter()
            .filter(|t| !single_decoy_ids.contains(t.compound_ref()))
            .cloned()
            .collect();
        exp.set_compounds(compounds);
        exp.set_transitions(transitions);
        Ok(())
    }

    /// Append unique ion signature (UIS) identification transitions for
    /// targets and decoys (source `uisTransitions`).
    ///
    /// 1. For every peptide, all alternative peptidoforms are enumerated
    ///    ([`theoretical_peptidoforms`](Self::theoretical_peptidoforms));
    ///    peptides with more than `max_num_alternative_localizations` are
    ///    skipped. Their theoretical ion series (and, with
    ///    `enable_ms2_precursors`, the precursor) are indexed per swath and
    ///    unmodified sequence.
    /// 2. Every distinct theoretical transition of a peptide that maps to at
    ///    least one peptidoform within `mz_threshold` becomes a non-detecting,
    ///    identifying, non-quantifying transition named
    ///    `<index>_UIS_{<peptidoforms>}_<precursor m/z>_<product m/z>_<RT>_<annotation>`,
    ///    with meta values `insilico_transition` (`"true"`) and `Peptidoforms`.
    /// 3. Unless `disable_decoy_transitions`, a random decoy sequence (seeded
    ///    by `shuffle_seed`) per unmodified target sequence receives the
    ///    target's modified residues, its peptidoforms are indexed the same
    ///    way, and decoy transitions named `<index>_UISDECOY_...` are appended,
    ///    except those overlapping a target peptidoform's ion.
    ///
    /// Only the transition list of `exp` changes; existing transitions stay
    /// first.
    ///
    /// # Errors
    ///
    /// Whatever peptide conversion, peptidoform enumeration and the ion
    /// series return; [`Error::InvalidValue`] when a UIS peptide has no
    /// retention time (the source's `getRetentionTime` throws), and when more
    /// than [`MAX_UIS_TRANSITIONS`] transitions would be appended. On error
    /// `exp` is unchanged.
    pub fn uis_transitions(exp: &mut TargetedExperiment, options: &UisOptions) -> Result<()> {
        let mut transitions = exp.transitions().to_vec();
        let mut target_sequences = SequenceMap::new();
        let mut target_ions = IonMap::new();
        let mut target_peptides = PeptideMap::new();
        generate_target_in_silico_map(
            exp,
            options,
            &mut target_sequences,
            &mut target_ions,
            &mut target_peptides,
        )?;
        sort_ion_map(&mut target_ions);
        sort_and_dedup_peptide_map(&mut target_peptides);
        generate_target_assays(
            exp,
            &mut transitions,
            options,
            &target_peptides,
            &target_ions,
        )?;

        if !options.disable_decoy_transitions {
            let seed = options.shuffle_seed.unwrap_or_else(time_seed);
            let decoy_sequences = generate_decoy_sequences(&target_sequences, seed)?;
            let mut target_decoy: BTreeMap<String, Peptide> = BTreeMap::new();
            let mut decoy_ions = IonMap::new();
            let mut decoy_peptides = PeptideMap::new();
            generate_decoy_in_silico_map(
                exp,
                options,
                &mut target_decoy,
                &target_peptides,
                &decoy_sequences,
                &mut decoy_ions,
                &mut decoy_peptides,
            )?;
            sort_ion_map(&mut decoy_ions);
            sort_and_dedup_peptide_map(&mut decoy_peptides);
            generate_decoy_assays(
                exp,
                &mut transitions,
                options,
                &decoy_peptides,
                &target_decoy,
                &decoy_ions,
                &target_ions,
            )?;
        }
        exp.set_transitions(transitions);
        Ok(())
    }
}

/// The modification name the source's `getModificationName` returns: the
/// record's short name, or a mass tag's full id.
fn modification_name(modification: &SequenceModification) -> String {
    modification.name().to_owned()
}

/// Group items by a string key in ascending key order, keeping input order
/// within a group (the source's `std::map<std::string, std::vector<...>>`).
fn group_by<'a, T>(
    items: &'a [T],
    key: impl Fn(&'a T) -> &'a str,
) -> BTreeMap<&'a str, Vec<&'a T>> {
    let mut groups: BTreeMap<&str, Vec<&T>> = BTreeMap::new();
    for item in items {
        groups.entry(key(item)).or_default().push(item);
    }
    groups
}

/// The bit patterns of the `max` largest library intensities of `group`
/// (source: sort with `std::greater`, resize, `std::set<double>`). Exact
/// equality, as the source's set lookup; `-0.0` and `0.0` are equal there and
/// are unified here.
fn top_intensities(group: &[&ReactionMonitoringTransition], max: usize) -> HashSet<u64> {
    let mut intensities: Vec<f64> = group.iter().map(|t| t.library_intensity()).collect();
    intensities.sort_by(|a, b| b.total_cmp(a));
    intensities.truncate(max);
    let mut set: HashSet<u64> = intensities.iter().map(|x| x.to_bits()).collect();
    if set.contains(&0.0_f64.to_bits()) || set.contains(&(-0.0_f64).to_bits()) {
        set.insert(0.0_f64.to_bits());
        set.insert((-0.0_f64).to_bits());
    }
    set
}

/// `std::set_difference` of two sorted multisets.
fn multiset_difference(left: &[String], right: &[String]) -> Vec<String> {
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < left.len() {
        if j >= right.len() || left[i] < right[j] {
            result.push(left[i].clone());
            i += 1;
        } else if right[j] < left[i] {
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    result
}

/// `std::prev_permutation` on a byte sequence; `false` (and the sequence
/// reset to ascending order) after the first permutation.
fn prev_permutation(values: &mut [u8]) -> bool {
    if values.len() < 2 {
        return false;
    }
    let mut i = values.len() - 1;
    while i > 0 && values[i - 1] <= values[i] {
        i -= 1;
    }
    if i == 0 {
        values.reverse();
        return false;
    }
    let mut j = values.len() - 1;
    while values[j] >= values[i - 1] {
        j -= 1;
    }
    values.swap(i - 1, j);
    values[i..].reverse();
    true
}

/// Whether `n choose k` exceeds `limit`.
fn binomial_exceeds(n: usize, k: usize, limit: usize) -> bool {
    let k = k.min(n - k);
    let mut value: u128 = 1;
    for i in 0..k {
        value = value * (n - i) as u128 / (i + 1) as u128;
        if value > limit as u128 {
            return true;
        }
    }
    false
}

fn swath_key(swath: i32) -> SwathKey {
    source_size(swath)
}

fn sort_ion_map(map: &mut IonMap) {
    for by_sequence in map.values_mut() {
        for ions in by_sequence.values_mut() {
            ions.sort_by(|a, b| a.0.total_cmp(&b.0));
        }
    }
}

fn sort_and_dedup_peptide_map(map: &mut PeptideMap) {
    for ions in map.values_mut() {
        ions.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
        ions.dedup();
    }
}

fn ion_series_for_uis(options: &UisOptions) -> IonSeriesOptions {
    let mut series = options.ion_series.clone();
    // Source: getIonSeries is called without round_decPow in the IPF maps.
    series.round_dec_pow = -4;
    series
}

fn peptide_sequence(peptide: &Peptide, options: &UisOptions) -> Result<AASequence> {
    Ok(get_aa_sequence_with(
        peptide,
        ModificationsDB::global(),
        options.ion_series.unresolved_modifications,
    )?
    .sequence)
}

fn generate_target_in_silico_map(
    exp: &TargetedExperiment,
    options: &UisOptions,
    sequences: &mut SequenceMap,
    ions: &mut IonMap,
    peptides: &mut PeptideMap,
) -> Result<()> {
    let series_options = ion_series_for_uis(options);
    let round = options.ion_series.round_dec_pow;
    for peptide in exp.peptides() {
        let sequence = peptide_sequence(peptide, options)?;
        let charge = peptide.charge_state().unwrap_or(1);
        let precursor_mz = typed_mz(&sequence, PeptideFragmentType::Full, charge)?;
        let swath = swath_key(MRMAssay::swath_index(&options.swathes, precursor_mz));
        let alternatives = MRMAssay::theoretical_peptidoforms(&sequence)?;
        if alternatives.len() > options.max_num_alternative_localizations {
            continue;
        }
        for alternative in &alternatives {
            let unmodified = alternative.as_str().to_owned();
            let modified = alternative.to_string();
            sequences
                .entry(swath)
                .or_default()
                .entry(unmodified.clone())
                .or_default()
                .insert(modified.clone());
            let series =
                MRMIonSeries::get_ion_series(alternative, source_size(charge), &series_options)?;
            index_peptidoform(
                ions,
                peptides,
                (swath, unmodified, modified),
                &peptide.id,
                &series,
                options.enable_ms2_precursors.then_some(precursor_mz),
                round,
            )?;
        }
    }
    Ok(())
}

/// Append one peptidoform's theoretical transitions to the ion index and the
/// peptide's transition list. As the source's `operator[]` calls, an entry is
/// created only when something is appended to it.
fn index_peptidoform(
    ions: &mut IonMap,
    peptides: &mut PeptideMap,
    (swath, unmodified, modified): (SwathKey, String, String),
    peptide_id: &str,
    series: &crate::analysis::openswath::mrm_ion_series::IonSeries,
    precursor_mz: Option<f64>,
    round: i32,
) -> Result<()> {
    let mut new_ions: FragmentSeqMap = Vec::new();
    let mut new_series: AssayIonSeries = Vec::new();
    if let Some(precursor_mz) = precursor_mz {
        let mz = round_decimal(precursor_mz, round)?;
        new_ions.push((mz, modified.clone()));
        new_series.push(("MS2_Precursor_i0".into(), mz));
    }
    for (annotation, &mz) in series {
        let mz = round_decimal(mz, round)?;
        new_ions.push((mz, modified.clone()));
        new_series.push((annotation.clone(), mz));
    }
    if !new_ions.is_empty() {
        ions.entry(swath)
            .or_default()
            .entry(unmodified)
            .or_default()
            .extend(new_ions);
        peptides
            .entry(peptide_id.to_owned())
            .or_default()
            .extend(new_series);
    }
    Ok(())
}

fn generate_decoy_sequences(
    target_sequences: &SequenceMap,
    seed: u32,
) -> Result<BTreeMap<String, String>> {
    let mut rng = BoostUniformInt::new(seed);
    let mut decoys: BTreeMap<String, String> = BTreeMap::new();
    for by_sequence in target_sequences.values() {
        for (unmodified, peptidoforms) in by_sequence {
            let mut decoy = match decoys.get(unmodified) {
                Some(existing) if !existing.is_empty() => existing.clone(),
                _ => MRMAssay::random_sequence(unmodified.len(), &mut rng),
            };
            decoys.entry(unmodified.clone()).or_default();
            for peptidoform in peptidoforms {
                let sequence = AASequence::parse(peptidoform)?;
                let residues = sequence.as_str().as_bytes();
                let mut bytes = decoy.into_bytes();
                let last = bytes.len().wrapping_sub(1);
                let mut copy = |i: usize| -> Result<()> {
                    match (bytes.get_mut(i), residues.get(i)) {
                        (Some(slot), Some(&residue)) => {
                            *slot = residue;
                            Ok(())
                        }
                        _ => Err(Error::InvalidValue(
                            "decoy and target sequence lengths differ".into(),
                        )),
                    }
                };
                if sequence.n_terminal_modification().is_some() {
                    copy(0)?;
                }
                if sequence.c_terminal_modification().is_some() {
                    copy(last)?;
                }
                for i in 0..sequence.len() {
                    if sequence.residue_modification(i)?.is_some() {
                        copy(i)?;
                    }
                }
                decoy = String::from_utf8_lossy(&bytes).into_owned();
                decoys.insert(unmodified.clone(), decoy.clone());
            }
        }
    }
    Ok(decoys)
}

fn generate_decoy_in_silico_map(
    exp: &TargetedExperiment,
    options: &UisOptions,
    target_decoy: &mut BTreeMap<String, Peptide>,
    target_peptides: &PeptideMap,
    decoy_sequences: &BTreeMap<String, String>,
    ions: &mut IonMap,
    peptides: &mut PeptideMap,
) -> Result<()> {
    let series_options = ion_series_for_uis(options);
    let round = options.ion_series.round_dec_pow;
    for peptide in exp.peptides() {
        let charge = peptide.charge_state().unwrap_or(1);
        if !target_peptides.contains_key(&peptide.id) {
            continue;
        }
        let sequence = peptide_sequence(peptide, options)?;
        let precursor_mz = typed_mz(&sequence, PeptideFragmentType::Full, charge)?;
        let swath = swath_key(MRMAssay::swath_index(&options.swathes, precursor_mz));
        let mut decoy_peptide = peptide.clone();
        decoy_peptide.sequence = decoy_sequences
            .get(&peptide.sequence)
            .cloned()
            .unwrap_or_default();
        target_decoy.insert(peptide.id.clone(), decoy_peptide.clone());
        let decoy_sequence = peptide_sequence(&decoy_peptide, options)?;
        let alternatives = MRMAssay::theoretical_peptidoforms_decoy(&sequence, &decoy_sequence)?;
        for alternative in &alternatives {
            let unmodified = alternative.as_str().to_owned();
            let modified = alternative.to_string();
            let series =
                MRMIonSeries::get_ion_series(alternative, source_size(charge), &series_options)?;
            // Source keys precursor entries by `peptide.id` and fragments by
            // `decoy_peptide.id`, a copy of the same id.
            index_peptidoform(
                ions,
                peptides,
                (swath, unmodified, modified),
                &decoy_peptide.id,
                &series,
                options.enable_ms2_precursors.then_some(precursor_mz),
                round,
            )?;
        }
    }
    Ok(())
}

fn ions_for<'a>(map: &'a IonMap, swath: SwathKey, unmodified: &str) -> Result<&'a FragmentSeqMap> {
    map.get(&swath)
        .and_then(|by_sequence| by_sequence.get(unmodified))
        .ok_or_else(|| {
            Error::MissingInformation(format!(
                "no in silico ions for {unmodified} in swath {swath}"
            ))
        })
}

fn push_uis(
    transitions: &mut Vec<ReactionMonitoringTransition>,
    tr: ReactionMonitoringTransition,
) -> Result<()> {
    if transitions.len() >= MAX_UIS_TRANSITIONS {
        return Err(Error::InvalidValue(format!(
            "UIS generation exceeds {MAX_UIS_TRANSITIONS} transitions"
        )));
    }
    transitions.push(tr);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn uis_transition(
    kind: &str,
    index: usize,
    isoforms: &[String],
    precursor_mz: f64,
    product_mz: f64,
    peptide: &Peptide,
    annotation: &str,
    decoy: bool,
    formula_loss: &crate::analysis::openswath::mrm_ion_series::FormulaLossParsing,
) -> Result<ReactionMonitoringTransition> {
    let mut trn = ReactionMonitoringTransition::new();
    if decoy {
        trn.set_decoy_transition_type(DecoyTransitionType::Decoy);
    }
    trn.set_detecting_transition(false);
    trn.cv_terms
        .metadata
        .insert("insilico_transition".into(), MetaValue::from("true"));
    trn.set_precursor_mz(precursor_mz);
    trn.set_product_mz(product_mz);
    trn.set_peptide_ref(peptide.id.clone());
    MRMIonSeries::annotate_transition_cv_with(&mut trn, annotation, formula_loss)?;
    trn.set_identifying_transition(true);
    trn.set_quantifying_transition(false);
    let joined = isoforms.join("|");
    let identifier = format!(
        "{index}_{kind}_{{{joined}}}_{}_{}_{}_{annotation}",
        format_float(trn.precursor_mz(), true),
        format_float(trn.product_mz(), true),
        format_float(peptide.retention_time()?, true),
    );
    trn.set_name(identifier.clone());
    trn.set_native_id(identifier);
    trn.cv_terms
        .metadata
        .insert("Peptidoforms".into(), MetaValue::from(joined));
    Ok(trn)
}

fn generate_target_assays(
    exp: &TargetedExperiment,
    transitions: &mut Vec<ReactionMonitoringTransition>,
    options: &UisOptions,
    target_peptides: &PeptideMap,
    target_ions: &IonMap,
) -> Result<()> {
    let round = options.ion_series.round_dec_pow;
    let mut index = 0usize;
    for (peptide_id, series) in target_peptides {
        let peptide = exp.peptide_by_ref(peptide_id)?;
        let charge = peptide.charge_state().unwrap_or(1);
        let sequence = peptide_sequence(peptide, options)?;
        let precursor_mz = typed_mz(&sequence, PeptideFragmentType::Full, charge)?;
        let swath = swath_key(MRMAssay::swath_index(&options.swathes, precursor_mz));
        let ions = ions_for(target_ions, swath, sequence.as_str())?;
        for (annotation, mz) in series {
            let isoforms = MRMAssay::matching_peptidoforms(*mz, ions, options.mz_threshold);
            if !isoforms.is_empty() {
                let trn = uis_transition(
                    "UIS",
                    index,
                    &isoforms,
                    round_decimal(precursor_mz, round)?,
                    *mz,
                    peptide,
                    annotation,
                    false,
                    &options.ion_series.formula_loss,
                )?;
                push_uis(transitions, trn)?;
            }
            index += 1;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn generate_decoy_assays(
    exp: &TargetedExperiment,
    transitions: &mut Vec<ReactionMonitoringTransition>,
    options: &UisOptions,
    decoy_peptides: &PeptideMap,
    target_decoy: &BTreeMap<String, Peptide>,
    decoy_ions: &IonMap,
    target_ions: &IonMap,
) -> Result<()> {
    let round = options.ion_series.round_dec_pow;
    let mut index = 0usize;
    for (peptide_id, series) in decoy_peptides {
        let target = exp.peptide_by_ref(peptide_id)?;
        let charge = target.charge_state().unwrap_or(1);
        let target_sequence = peptide_sequence(target, options)?;
        let precursor_mz = typed_mz(&target_sequence, PeptideFragmentType::Full, charge)?;
        let swath = swath_key(MRMAssay::swath_index(&options.swathes, precursor_mz));
        let decoy_peptide = target_decoy.get(peptide_id).cloned().unwrap_or_default();
        let decoy_sequence = peptide_sequence(&decoy_peptide, options)?;
        let decoy_ion_list = ions_for(decoy_ions, swath, decoy_sequence.as_str())?;
        let target_ion_list = ions_for(target_ions, swath, target_sequence.as_str())?;
        for (annotation, mz) in series {
            let isoforms =
                MRMAssay::matching_peptidoforms(*mz, decoy_ion_list, options.mz_threshold);
            if !isoforms.is_empty() {
                let trn = uis_transition(
                    "UISDECOY",
                    index,
                    &isoforms,
                    round_decimal(precursor_mz, round)?,
                    *mz,
                    &decoy_peptide,
                    annotation,
                    true,
                    &options.ion_series.formula_loss,
                )?;
                let overlap =
                    MRMAssay::matching_peptidoforms(*mz, target_ion_list, options.mz_threshold);
                if !overlap.is_empty() {
                    // Source: `continue` skips the index increment too.
                    continue;
                }
                push_uis(transitions, trn)?;
            }
            index += 1;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prev_permutation_matches_the_standard_library_contract() {
        let mut v = vec![1u8, 1, 0];
        assert!(prev_permutation(&mut v));
        assert_eq!(v, vec![1, 0, 1]);
        assert!(prev_permutation(&mut v));
        assert_eq!(v, vec![0, 1, 1]);
        assert!(!prev_permutation(&mut v));
        assert_eq!(v, vec![1, 1, 0]);
    }

    #[test]
    fn multiset_difference_keeps_surplus_copies() {
        let owned = |v: &[&str]| v.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>();
        assert_eq!(
            multiset_difference(&owned(&["a", "a", "b", "c"]), &owned(&["a", "c"])),
            owned(&["a", "b"])
        );
    }

    #[test]
    fn binomial_limit() {
        assert!(!binomial_exceeds(16, 5, 4368));
        assert!(binomial_exceeds(16, 5, 4367));
        assert!(!binomial_exceeds(3, 0, 1));
    }
}
