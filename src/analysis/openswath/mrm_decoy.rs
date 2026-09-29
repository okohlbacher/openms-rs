// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Decoy generation for targeted assay libraries.
//!
//! Ports `ANALYSIS/OPENSWATH/MRMDecoy.h` (with its `.cpp`), the algorithm of
//! the `OpenSwathDecoyGenerator` TOPP tool. [`MRMDecoy::generate_decoys`]
//! builds a decoy [`TargetedExperiment`] from a target one.
//!
//! The simplest methods are *reverse* and *pseudo-reverse*, which reverse the
//! peptide sequence completely or leave the last (tryptic) residue in place.
//! *Shuffle* uses an algorithm similar to the one of Lam et al. (2010),
//! "Artificial decoy spectral libraries for false discovery rate estimation in
//! spectral library searching in proteomics", J. Proteome Res. 9, 605-610: it
//! shuffles the sequence and carries the fragment intensities over, which
//! needs the fragment ions annotated first.
//!
//! First the decoy method is applied to every target peptide. Then, per
//! peptide, each target fragment ion is matched to its most likely origin (ion
//! series, ordinal and optionally one of the neutral losses of
//! [`MRMIonSeries`]), and the equivalent ion of the decoy peptide is created
//! with the same library intensity: a target `b5` with intensity 200 becomes a
//! decoy `b5` with intensity 200. The decoy product m/z is the theoretical
//! value.
//!
//! # Not ported
//!
//! The source's `LightTargetedExperiment` overloads (`generateDecoysLight`,
//! `reversePeptideLight`, `shufflePeptideLight`, `pseudoreversePeptideLight_`,
//! `switchKRLight`, `hasCNterminalModsLight_`) operate on
//! `OPENSWATHALGO/DATAACCESS/TransitionExperiment.h`, which is not ported yet.
//! The heavy methods here delegate to the same string-level algorithms the
//! source's light methods implement, so the light overloads reduce to thin
//! adapters once that header lands. See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/MRM_ASSAY_SUPPORT.md).
//!
//! The source contains no `#pragma omp`; the port is serial.
//!
//! [`MRMDecoy::generate_decoys`]: crate::analysis::openswath::mrm_decoy::MRMDecoy::generate_decoys
//! [`MRMIonSeries`]: crate::analysis::openswath::mrm_ion_series::MRMIonSeries
//! [`TargetedExperiment`]: crate::analysis::openswath::targeted_experiment::TargetedExperiment

use crate::analysis::openswath::mrm_ion_series::{
    IonSeriesOptions, MRMIonSeries, source_size, typed_mz,
};
use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{
    Peptide, Protein, get_aa_sequence_with,
};
use crate::chemistry::decoy_generator::DecoyRandom;
use crate::chemistry::{ModificationsDB, PeptideFragmentType};
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, Result};
use rand_mt::Mt;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::{SystemTime, UNIX_EPOCH};

/// The source `DefaultParamHandler` name.
pub const MRM_DECOY_NAME: &str = "MRMDecoy";

/// Residue positions, as the source `IndexType` (`std::vector<size_t>`).
pub type IndexType = Vec<usize>;

/// Native ceiling on `max_attempts * sequence length` for one shuffle: 2^28.
pub const MAX_SHUFFLE_WORK: usize = 1 << 28;

/// Native ceiling on the number of transitions [`MRMDecoy::generate_decoys`]
/// may produce: 2^26.
pub const MAX_DECOY_TRANSITIONS: usize = 1 << 26;

/// The 17 residues the source draws random replacements from, in source order.
pub(crate) const RANDOM_RESIDUES: [u8; 17] = *b"ANDCEQGHILMFSTWYV";

/// `boost::variate_generator<boost::mt19937&, boost::uniform_int<>>`, the
/// source's random source for shuffling and random decoy sequences.
///
/// Raw words come from `rand_mt::Mt`, the reference MT19937, seeded as
/// `boost::mt19937(seed)`. `boost::uniform_int<>` defaults to the range
/// `[0, 9]`, so [`draw`](Self::draw) (the source's `pseudoRNG()`) yields at
/// most 9; [`below`](Self::below) (`pseudoRNG(n)`) yields `[0, n - 1]`. Both
/// map a word with Boost's bucket division and rejection, not a modulo.
#[derive(Clone, Debug)]
pub(crate) struct BoostUniformInt {
    engine: Mt,
}

impl BoostUniformInt {
    /// The generator in the state of `boost::mt19937(seed)`.
    pub(crate) fn new(seed: u32) -> Self {
        Self {
            engine: Mt::new(seed),
        }
    }

    /// `pseudoRNG()`: uniform in `[0, 9]`.
    pub(crate) fn draw(&mut self) -> u32 {
        self.bounded(9)
    }

    /// `pseudoRNG(n)`: uniform in `[0, n - 1]`; `n` of 0 or 1 draws nothing.
    pub(crate) fn below(&mut self, n: u32) -> u32 {
        self.bounded(n.saturating_sub(1))
    }

    // Boost `generate_uniform_int` for a 32-bit engine: bucket size
    // floor(2^32 / (high + 1)), words past the last full bucket rejected.
    fn bounded(&mut self, high: u32) -> u32 {
        if high == 0 {
            return 0;
        }
        if high == u32::MAX {
            return self.engine.next_u32();
        }
        let bucket = ((1_u64 << 32) / (u64::from(high) + 1)) as u32;
        loop {
            let result = self.engine.next_u32() / bucket;
            if result <= high {
                return result;
            }
        }
    }
}

/// The seed the source derives from `time(nullptr)` when a seed of `-1` is
/// passed: the low 32 bits of the current Unix time in seconds.
pub(crate) fn time_seed() -> u32 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() as u32)
}

/// The decoy method of [`MRMDecoy::generate_decoys`] (the source `method`
/// string).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DecoyMethod {
    /// `"shuffle"`: shuffle the non-fixed residues.
    #[default]
    Shuffle,
    /// `"reverse"`: reverse the whole sequence.
    Reverse,
    /// `"pseudo-reverse"`: reverse all but the C-terminal residue.
    PseudoReverse,
    /// `"shift"`: keep the sequence and shift the product m/z by
    /// [`DecoyOptions::product_mz_shift`]. As the source, the unchanged decoy
    /// precursor then equals its target and every such decoy is excluded (see
    /// the support document).
    Shift,
}

impl DecoyMethod {
    /// Parse the source method string.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] for any other string. The source silently
    /// treats an unknown method like `"shift"` without the shift; the
    /// `OpenSwathDecoyGenerator` tool restricts the string to these four.
    pub fn from_source_str(method: &str) -> Result<Self> {
        match method {
            "shuffle" => Ok(Self::Shuffle),
            "reverse" => Ok(Self::Reverse),
            "pseudo-reverse" => Ok(Self::PseudoReverse),
            "shift" => Ok(Self::Shift),
            _ => Err(Error::InvalidValue(format!(
                "unknown decoy method '{method}' (valid: shuffle, reverse, pseudo-reverse, shift)"
            ))),
        }
    }
}

/// The arguments of the source `generateDecoys` after the two experiments.
#[derive(Clone, Debug, PartialEq)]
pub struct DecoyOptions {
    /// Decoy method.
    pub method: DecoyMethod,
    /// Fraction of target peptides to generate decoys for, at most 1. Below 1
    /// a deterministic shuffle (the source `Math::RandomShuffler` with its
    /// default seed) selects them, cycling when the fraction needs more.
    pub aim_decoy_fraction: f64,
    /// Switch the terminal K/R (K to R and R to K; another residue becomes a
    /// deterministic pseudo-random one), so the decoy precursor and y ions
    /// differ in mass. See Bruderer et al., Mol Cell Proteomics 2017,
    /// 10.1074/mcp.RA117.000314 (supplemental material).
    pub switch_kr: bool,
    /// Prefix of decoy protein, peptide and transition identifiers.
    pub decoy_tag: String,
    /// Maximum shuffle attempts per peptide.
    pub max_attempts: usize,
    /// Maximum sequence identity of a shuffled decoy to its target.
    pub identity_threshold: f64,
    /// Added to every decoy precursor m/z, so TOPPView can tell decoy and
    /// target precursors apart.
    pub precursor_mz_shift: f64,
    /// Added to every decoy product m/z with [`DecoyMethod::Shift`].
    pub product_mz_shift: f64,
    /// Threshold for matching a target product m/z to its theoretical ion.
    pub product_mz_threshold: f64,
    /// Fragment types, charges, losses and rounding of the ion series.
    pub ion_series: IonSeriesOptions,
    /// Shuffle seed. `None` seeds every peptide's shuffle from the current
    /// time in seconds, as the source always does (so peptides shuffled within
    /// the same second share a seed); `Some` makes the output reproducible.
    pub shuffle_seed: Option<u32>,
}

impl DecoyOptions {
    /// The `OpenSwathDecoyGenerator` defaults: shuffle, fraction 1, no K/R
    /// switch, tag `DECOY_`, 30 attempts, identity 0.7, precursor shift 0,
    /// product shift 20, product threshold 0.025, and `ion_series` as given.
    pub fn new(ion_series: IonSeriesOptions) -> Self {
        Self {
            method: DecoyMethod::Shuffle,
            aim_decoy_fraction: 1.0,
            switch_kr: false,
            decoy_tag: "DECOY_".into(),
            max_attempts: 30,
            identity_threshold: 0.7,
            precursor_mz_shift: 0.0,
            product_mz_shift: 20.0,
            product_mz_threshold: 0.025,
            ion_series,
            shuffle_seed: None,
        }
    }
}

/// Generates decoy peptides and transitions from a [`TargetedExperiment`].
///
/// Parameters (source defaults): `non_shuffle_pattern` `"KRP"`, the residues
/// that keep their position when shuffling (lysine, arginine and proline);
/// `keepPeptideNTerm` and `keepPeptideCTerm` `"true"` (both advanced), whether
/// the peptide termini stay in place when shuffling. Reversal ignores all
/// three, as the source's `reversePeptide_` and `pseudoreversePeptide_` do.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MRMDecoy {
    /// Residues not shuffled (`non_shuffle_pattern`).
    pub keep_const_pattern: String,
    /// Keep the N-terminal residue in place when shuffling
    /// (`keepPeptideNTerm`).
    pub keep_n: bool,
    /// Keep the C-terminal residue in place when shuffling
    /// (`keepPeptideCTerm`).
    pub keep_c: bool,
}

impl Default for MRMDecoy {
    fn default() -> Self {
        Self {
            keep_const_pattern: "KRP".into(),
            keep_n: true,
            keep_c: true,
        }
    }
}

fn entry(name: &str, value: ParamValue, description: &str) -> ParamEntry {
    ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    }
}

fn flag_entry(name: &str, value: bool, description: &str) -> ParamEntry {
    let mut e = entry(
        name,
        ParamValue::String(if value { "true" } else { "false" }.into()),
        description,
    );
    e.valid_strings = vec!["true".into(), "false".into()];
    e.tags.insert("advanced".into());
    e
}

impl MRMDecoy {
    /// A decoy generator with the source defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// The source defaults as a [`Param`] (`getDefaults`).
    ///
    /// # Errors
    ///
    /// Propagates a [`Param`] construction failure, which does not occur for
    /// these three fixed entries.
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The current parameters in the layout of [`defaults`](Self::defaults)
    /// (`getParameters`).
    ///
    /// # Errors
    ///
    /// As [`defaults`](Self::defaults).
    pub fn to_param(&self) -> Result<Param> {
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![
                entry(
                    "non_shuffle_pattern",
                    ParamValue::String(self.keep_const_pattern.clone()),
                    "Residues to not shuffle (keep at a constant position when shuffling). Default is 'KPR' to not shuffle lysine, arginine and proline.",
                ),
                flag_entry(
                    "keepPeptideNTerm",
                    self.keep_n,
                    "Whether to keep peptide N terminus constant when shuffling / reversing.",
                ),
                flag_entry(
                    "keepPeptideCTerm",
                    self.keep_c,
                    "Whether to keep peptide C terminus constant when shuffling / reversing.",
                ),
            ],
            nodes: Vec::new(),
        })
    }

    /// Build a generator from parameters, as the source `setParameters`
    /// followed by `updateMembers_`. Missing parameters take their defaults;
    /// unknown ones are ignored (see
    /// [`from_param_with_warnings`](Self::from_param_with_warnings)).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a parameter has the wrong type or
    /// violates its valid strings (the source `Exception::InvalidParameter`).
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`from_param`](Self::from_param), also returning the
    /// unknown-parameter warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`from_param`](Self::from_param).
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(MRM_DECOY_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, |param| {
            Ok(Self {
                keep_const_pattern: param.value("non_shuffle_pattern")?.as_str()?.to_owned(),
                keep_n: param.value("keepPeptideNTerm")?.to_bool()?,
                keep_c: param.value("keepPeptideCTerm")?.to_bool()?,
            })
        })
    }

    /// The positions of `sequence` that are not reversed or shuffled.
    ///
    /// A position is fixed when it is the first and `keep_n` is set, the last
    /// and `keep_c` is set, or its residue occurs in `keep_const_pattern`. As
    /// the source, a non-terminal position is listed once per matching pattern
    /// character, so a pattern with a repeated letter lists it repeatedly.
    pub fn find_fixed_residues(
        sequence: &str,
        keep_n: bool,
        keep_c: bool,
        keep_const_pattern: &str,
    ) -> IndexType {
        let bytes = sequence.as_bytes();
        let mut idx = Vec::new();
        for (i, &residue) in bytes.iter().enumerate() {
            if (keep_n && i == 0) || (keep_c && i + 1 == bytes.len()) {
                idx.push(i);
                continue;
            }
            for &fixed in keep_const_pattern.as_bytes() {
                if residue == fixed {
                    idx.push(i);
                }
            }
        }
        idx
    }

    /// The positions fixed by `keep_const_pattern` alone (source
    /// `findFixedResidues_`, which the SpectraST decoy generator inspired).
    pub fn fixed_residues(&self, sequence: &str) -> IndexType {
        Self::find_fixed_residues(sequence, false, false, &self.keep_const_pattern)
    }

    /// The positions fixed by `keep_const_pattern` and the configured termini
    /// (source `findFixedAndTermResidues_`).
    pub fn fixed_and_terminal_residues(&self, sequence: &str) -> IndexType {
        Self::find_fixed_residues(sequence, self.keep_n, self.keep_c, &self.keep_const_pattern)
    }

    /// Relative number of positions at which `sequence` and `decoy` hold the
    /// same residue (source `AASequenceIdentity`).
    ///
    /// The source computes a `double` and returns it as `float`; so does this.
    /// Two empty sequences give NaN (`0 / 0`), as the source.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when the lengths differ (the source's debug-only
    /// precondition; its release build reads past the shorter one).
    pub fn aa_sequence_identity(sequence: &str, decoy: &str) -> Result<f32> {
        if sequence.len() != decoy.len() {
            return Err(Error::InvalidValue(
                "cannot compare two sequences of unequal length".into(),
            ));
        }
        let matches = sequence
            .bytes()
            .zip(decoy.bytes())
            .filter(|(a, b)| a == b)
            .count();
        Ok((matches as f64 / sequence.len() as f64) as f32)
    }

    /// Shuffle `peptide`'s sequence, carrying its modifications, until the
    /// identity to the (possibly mutated) original is at most
    /// `identity_threshold` or `max_attempts` shuffles were made.
    ///
    /// Each attempt shuffles the positions not fixed by
    /// [`fixed_and_terminal_residues`](Self::fixed_and_terminal_residues)
    /// with a forward Fisher-Yates over `boost::mt19937(seed)`, then relocates
    /// every residue modification with its residue (terminal modifications at
    /// `-1` and the sequence length stay). After attempts 9, 19, 29, ... one
    /// random unmodified, non-terminal residue is replaced by one of the first
    /// ten of `ANDCEQGHILMFSTWYV` (the source draws `pseudoRNG() % 17` from a
    /// `[0, 9]` distribution) and the mutation is kept for later attempts.
    ///
    /// The source's `seed = -1` (seed from the current time) is
    /// [`DecoyOptions::shuffle_seed`] `None` in
    /// [`generate_decoys`](Self::generate_decoys); here the seed is explicit.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when `max_attempts * sequence length` exceeds
    /// [`MAX_SHUFFLE_WORK`] (native), or when a repeated letter in
    /// `keep_const_pattern` makes the source's index erasure run out of range.
    pub fn shuffle_peptide(
        &self,
        peptide: &Peptide,
        identity_threshold: f64,
        seed: u32,
        max_attempts: usize,
    ) -> Result<Peptide> {
        let locations: Vec<i32> = peptide.mods.iter().map(|m| m.location).collect();
        let (sequence, locations) = self.shuffle_sequence(
            &peptide.sequence,
            &locations,
            identity_threshold,
            seed,
            max_attempts,
        )?;
        let mut shuffled = peptide.clone();
        shuffled.sequence = sequence;
        for (modification, location) in shuffled.mods.iter_mut().zip(locations) {
            modification.location = location;
        }
        Ok(shuffled)
    }

    /// The string-level shuffle of the source `shufflePeptideLight`, on a
    /// sequence and its modification locations.
    fn shuffle_sequence(
        &self,
        sequence: &str,
        locations: &[i32],
        identity_threshold: f64,
        seed: u32,
        max_attempts: usize,
    ) -> Result<(String, Vec<i32>)> {
        if max_attempts
            .checked_mul(sequence.len())
            .is_none_or(|work| work > MAX_SHUFFLE_WORK)
        {
            return Err(Error::InvalidValue(format!(
                "shuffle work max_attempts * length exceeds {MAX_SHUFFLE_WORK}"
            )));
        }
        let mut shuffled: Vec<u8> = sequence.as_bytes().to_vec();
        let mut shuffled_locations = locations.to_vec();
        let mut peptide_seq: Vec<u8> = shuffled.clone();
        let mut rng = BoostUniformInt::new(seed);
        let identity = |a: &[u8], b: &[u8]| -> f64 {
            let matches = a.iter().zip(b).filter(|(x, y)| x == y).count();
            f64::from((matches as f64 / a.len() as f64) as f32)
        };
        let mut attempts = 0usize;
        while identity(&peptide_seq, &shuffled) > identity_threshold && attempts < max_attempts {
            let text = String::from_utf8_lossy(&peptide_seq).into_owned();
            let idx = self.fixed_and_terminal_residues(&text);
            shuffled.clone_from(&peptide_seq);
            shuffled_locations = locations.to_vec();

            let mut peptide_index: Vec<usize> = (0..peptide_seq.len()).collect();
            remove_fixed(&mut peptide_index, &idx)?;
            for i in 1..peptide_index.len() {
                let n = u32::try_from(i + 1)
                    .map_err(|_| Error::InvalidValue("sequence too long to shuffle".into()))?;
                let j = rng.below(n) as usize;
                peptide_index.swap(i, j);
            }
            reinsert_fixed(&mut peptide_index, &idx)?;
            for (i, &from) in peptide_index.iter().enumerate() {
                shuffled[i] = peptide_seq[from];
            }
            relocate(&mut shuffled_locations, &peptide_index, peptide_seq.len());

            attempts += 1;
            if attempts % 10 == 9 {
                // Source: pick the residue first, then the position.
                let residue = RANDOM_RESIDUES[(rng.draw() % 17) as usize];
                let len = shuffled.len();
                let mut trials = 0usize;
                while trials < len {
                    let position = (rng.draw() as usize) % len;
                    trials += 1;
                    let modified = shuffled_locations
                        .iter()
                        .any(|&l| usize::try_from(l).is_ok_and(|l| l == position));
                    if !(modified || position == 0 || position == len - 1) {
                        shuffled[position] = residue;
                        break;
                    }
                }
                peptide_seq.clone_from(&shuffled);
            }
        }
        Ok((
            String::from_utf8_lossy(&shuffled).into_owned(),
            shuffled_locations,
        ))
    }

    /// Reverse `peptide`'s sequence, carrying its modifications.
    ///
    /// Positions fixed by [`find_fixed_residues`](Self::find_fixed_residues)
    /// with `keep_n`, `keep_c` and `const_pattern` keep their residue; the
    /// others are reversed among themselves. Residue modifications move with
    /// their residue; terminal modifications (`-1`, sequence length) stay.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a repeated letter in `const_pattern` makes
    /// the source's index erasure run out of range (undefined behaviour in the
    /// source).
    pub fn reverse_peptide(
        peptide: &Peptide,
        keep_n: bool,
        keep_c: bool,
        const_pattern: &str,
    ) -> Result<Peptide> {
        let idx = Self::find_fixed_residues(&peptide.sequence, keep_n, keep_c, const_pattern);
        let bytes = peptide.sequence.as_bytes();
        let mut peptide_index: Vec<usize> = (0..bytes.len()).collect();
        remove_fixed(&mut peptide_index, &idx)?;
        peptide_index.reverse();
        reinsert_fixed(&mut peptide_index, &idx)?;
        let reversed: Vec<u8> = peptide_index.iter().map(|&i| bytes[i]).collect();
        let mut locations: Vec<i32> = peptide.mods.iter().map(|m| m.location).collect();
        relocate(&mut locations, &peptide_index, bytes.len());
        let mut result = peptide.clone();
        result.sequence = String::from_utf8_lossy(&reversed).into_owned();
        for (modification, location) in result.mods.iter_mut().zip(locations) {
            modification.location = location;
        }
        Ok(result)
    }

    /// Pseudo-reverse: reverse all but the C-terminal residue, with no
    /// pattern (source `pseudoreversePeptide_`).
    ///
    /// # Errors
    ///
    /// As [`reverse_peptide`](Self::reverse_peptide), which cannot fail here.
    pub fn pseudoreverse_peptide(peptide: &Peptide) -> Result<Peptide> {
        Self::reverse_peptide(peptide, false, true, "")
    }

    /// Reverse the whole sequence, with no pattern (source `reversePeptide_`).
    ///
    /// # Errors
    ///
    /// As [`reverse_peptide`](Self::reverse_peptide), which cannot fail here.
    pub fn reverse_peptide_full(peptide: &Peptide) -> Result<Peptide> {
        Self::reverse_peptide(peptide, false, false, "")
    }

    /// Switch the final residue of a tryptic peptide: K becomes R and R
    /// becomes K.
    ///
    /// Any other final residue is replaced by one of `ANDCEQGHILMFSTWYV`
    /// chosen by the 64-bit FNV-1a hash of the sequence modulo 17, so the same
    /// peptide always gets the same replacement. The source hashes `char`
    /// values; bytes above 127 are sign-extended here as on x86-64, where the
    /// source's `char` is signed. An empty sequence is left alone.
    pub fn switch_kr(peptide: &mut Peptide) {
        let mut bytes = std::mem::take(&mut peptide.sequence).into_bytes();
        if let Some(last) = bytes.len().checked_sub(1) {
            bytes[last] = match bytes[last] {
                b'K' => b'R',
                b'R' => b'K',
                _ => {
                    let mut hash: u64 = 14_695_981_039_346_656_037;
                    for &byte in &bytes {
                        hash ^= i64::from(byte as i8) as u64;
                        hash = hash.wrapping_mul(1_099_511_628_211);
                    }
                    RANDOM_RESIDUES[(hash % 17) as usize]
                }
            };
        }
        peptide.sequence = String::from_utf8_lossy(&bytes).into_owned();
    }

    /// Whether `peptide` has an N-terminal (`-1`) or C-terminal (sequence
    /// length) modification, or, with `check_c_terminal_aa`, a modification of
    /// its last residue (source `hasCNterminalMods_`).
    pub fn has_cn_terminal_mods(peptide: &Peptide, check_c_terminal_aa: bool) -> bool {
        let len = i64::try_from(peptide.sequence.len()).unwrap_or(i64::MAX);
        peptide.mods.iter().any(|m| {
            let location = i64::from(m.location);
            location == -1 || location == len || (check_c_terminal_aa && location == len - 1)
        })
    }

    /// The peptide sequence with `(UniMod:<id>)` after each modified position,
    /// N-terminal ones first and C-terminal ones last (source
    /// `getModifiedPeptideSequence_`).
    pub fn modified_peptide_sequence(peptide: &Peptide) -> String {
        let bytes = peptide.sequence.as_bytes();
        let mut result = String::new();
        let append_mods = |result: &mut String, location: i64| {
            for m in &peptide.mods {
                if i64::from(m.location) == location {
                    result.push_str(&format!("(UniMod:{})", m.unimod_id));
                }
            }
        };
        append_mods(&mut result, -1);
        for (i, &residue) in bytes.iter().enumerate() {
            result.push(char::from(residue));
            append_mods(&mut result, i as i64);
        }
        append_mods(&mut result, bytes.len() as i64);
        result
    }

    /// Generate decoys for the peptides and transitions of `exp` into `dec`.
    ///
    /// Proteins, peptides and peptide group labels get `options.decoy_tag`
    /// prepended to their ids. Every selected peptide is converted with
    /// `options.method`; with [`DecoyMethod::PseudoReverse`] and
    /// [`DecoyMethod::Reverse`] a peptide with terminal modifications (for
    /// pseudo-reverse with `switch_kr`, also a modified last residue) is
    /// excluded, and with [`DecoyMethod::Shuffle`] and `switch_kr` a shuffled
    /// peptide with such modifications is excluded. With `switch_kr` the last
    /// residue is then switched ([`switch_kr`](Self::switch_kr)). A decoy
    /// whose modified sequence and charge equal a selected target or an
    /// earlier decoy is excluded.
    ///
    /// Transitions are visited grouped by peptide reference in ascending
    /// order. For every detecting, non-decoy transition of a peptide with a
    /// decoy, the target product m/z is annotated against the target ion
    /// series ([`MRMIonSeries::annotate_ion`] within
    /// `options.product_mz_threshold`) and the same ion is looked up in the
    /// decoy series. The decoy transition copies the target one with native id
    /// and peptide reference tagged, type decoy, the decoy precursor m/z
    /// (theoretical, plus `options.precursor_mz_shift`) and the decoy ion m/z
    /// (plus `options.product_mz_shift` for [`DecoyMethod::Shift`]). A
    /// transition whose ion has no decoy counterpart excludes its whole
    /// peptide. Finally `dec` receives the decoy transitions, peptides and
    /// proteins that survived the exclusions; its other content is kept.
    ///
    /// # Errors
    ///
    /// - [`Error::InvalidValue`] when `options.aim_decoy_fraction` exceeds 1
    ///   (the source `Exception::IllegalArgument`), is negative or NaN, or when
    ///   more than [`MAX_DECOY_TRANSITIONS`] transitions would result (native).
    /// - Whatever peptide conversion, the ion series, the reversal or the
    ///   shuffle return, and [`Error::MissingInformation`] when a transition's
    ///   peptide is missing.
    ///
    /// On error `dec` is unchanged.
    pub fn generate_decoys(
        &self,
        exp: &TargetedExperiment,
        dec: &mut TargetedExperiment,
        options: &DecoyOptions,
    ) -> Result<()> {
        let fraction = options.aim_decoy_fraction;
        if fraction > 1.0 {
            return Err(Error::InvalidValue(
                "Decoy fraction needs to be less than one (values larger than one currently not supported).".into(),
            ));
        }
        if !(fraction >= 0.0) {
            return Err(Error::InvalidValue(
                "decoy fraction must not be negative or NaN".into(),
            ));
        }
        let tag = options.decoy_tag.as_str();
        let proteins: Vec<Protein> = exp
            .proteins()
            .iter()
            .map(|protein| {
                let mut protein = protein.clone();
                protein.id = format!("{tag}{}", protein.id);
                protein
            })
            .collect();

        let n = exp.peptides().len();
        let mut selection: Vec<usize> = (0..n).collect();
        if fraction < 1.0 {
            let mut shuffler = DecoyRandom::seeded(5489);
            let mut draws = n.saturating_mul(64).max(64);
            shuffler.shuffle_slice(&mut selection, &mut draws)?;
            let item_list = selection;
            selection = Vec::new();
            let mut k = 0usize;
            while (selection.len() as f64) < fraction * n as f64 && !item_list.is_empty() {
                selection.push(item_list[k % item_list.len()]);
                k += 1;
            }
        }

        let policy = options.ion_series.unresolved_modifications;
        let mut all_sequences: HashMap<String, String> = HashMap::new();
        for &index in &selection {
            let peptide = &exp.peptides()[index];
            all_sequences.insert(precursor_key(peptide), peptide.id.clone());
        }

        let mut exclusion: HashSet<String> = HashSet::new();
        let mut peptides: Vec<Peptide> = Vec::with_capacity(selection.len());
        for &index in &selection {
            let mut peptide = exp.peptides()[index].clone();
            peptide.id = format!("{tag}{}", peptide.id);
            if !peptide.peptide_group_label().is_empty() {
                let label = format!("{tag}{}", peptide.peptide_group_label());
                peptide.set_peptide_group_label(label);
            }
            match options.method {
                DecoyMethod::PseudoReverse => {
                    if Self::has_cn_terminal_mods(&peptide, options.switch_kr) {
                        exclusion.insert(peptide.id.clone());
                    } else {
                        peptide = Self::pseudoreverse_peptide(&peptide)?;
                        if options.switch_kr {
                            Self::switch_kr(&mut peptide);
                        }
                    }
                }
                DecoyMethod::Reverse => {
                    if Self::has_cn_terminal_mods(&peptide, false) {
                        exclusion.insert(peptide.id.clone());
                    } else {
                        peptide = Self::reverse_peptide_full(&peptide)?;
                    }
                }
                DecoyMethod::Shuffle => {
                    let seed = options.shuffle_seed.unwrap_or_else(time_seed);
                    peptide = self.shuffle_peptide(
                        &peptide,
                        options.identity_threshold,
                        seed,
                        options.max_attempts,
                    )?;
                    if options.switch_kr && Self::has_cn_terminal_mods(&peptide, true) {
                        exclusion.insert(peptide.id.clone());
                    } else if options.switch_kr {
                        Self::switch_kr(&mut peptide);
                    }
                }
                DecoyMethod::Shift => {}
            }
            let key = precursor_key(&peptide);
            if all_sequences.contains_key(&key) {
                exclusion.insert(peptide.id.clone());
            } else {
                all_sequences.insert(key, peptide.id.clone());
            }
            for protein_ref in &mut peptide.protein_refs {
                *protein_ref = format!("{tag}{protein_ref}");
            }
            peptides.push(peptide);
        }
        // Source `dec.setPeptides(peptides)` then `dec.hasPeptide` /
        // `getPeptideByRef`: the last peptide with an id wins.
        let decoy_index: HashMap<&str, usize> = peptides
            .iter()
            .enumerate()
            .map(|(i, p)| (p.id.as_str(), i))
            .collect();

        let mut by_peptide: BTreeMap<&str, Vec<&ReactionMonitoringTransition>> = BTreeMap::new();
        for transition in exp.transitions() {
            by_peptide
                .entry(transition.peptide_ref())
                .or_default()
                .push(transition);
        }

        let mut transitions: Vec<ReactionMonitoringTransition> = Vec::new();
        for (peptide_ref, group) in &by_peptide {
            let decoy_ref = format!("{tag}{peptide_ref}");
            let Some(&decoy_position) = decoy_index.get(decoy_ref.as_str()) else {
                continue;
            };
            let target = exp.peptide_by_ref(peptide_ref)?;
            let decoy = &peptides[decoy_position];
            let db = ModificationsDB::global();
            let target_sequence = get_aa_sequence_with(target, db, policy)?.sequence;
            let decoy_sequence = get_aa_sequence_with(decoy, db, policy)?.sequence;
            let decoy_charge = decoy.charge_state().unwrap_or(1);
            let target_charge = target.charge_state().unwrap_or(1);
            let decoy_series = MRMIonSeries::get_ion_series(
                &decoy_sequence,
                source_size(decoy_charge),
                &options.ion_series,
            )?;
            let target_series = MRMIonSeries::get_ion_series(
                &target_sequence,
                source_size(target_charge),
                &options.ion_series,
            )?;
            let decoy_precursor_mz =
                typed_mz(&decoy_sequence, PeptideFragmentType::Full, decoy_charge)?
                    + options.precursor_mz_shift;

            for &tr in group {
                if !tr.is_detecting_transition()
                    || tr.decoy_transition_type() == DecoyTransitionType::Decoy
                {
                    continue;
                }
                let mut decoy_tr = tr.clone();
                decoy_tr.set_native_id(format!("{tag}{}", tr.native_id()));
                decoy_tr.set_decoy_transition_type(DecoyTransitionType::Decoy);
                decoy_tr.set_precursor_mz(decoy_precursor_mz);
                let (target_ion, _) = MRMIonSeries::annotate_ion(
                    &target_series,
                    tr.product_mz(),
                    options.product_mz_threshold,
                );
                let (_, decoy_mz) = MRMIonSeries::get_ion(&decoy_series, &target_ion);
                if options.method == DecoyMethod::Shift {
                    decoy_tr.set_product_mz(decoy_mz + options.product_mz_shift);
                } else {
                    decoy_tr.set_product_mz(decoy_mz);
                }
                decoy_tr.set_peptide_ref(format!("{tag}{}", tr.peptide_ref()));
                if decoy_mz > 0.0 {
                    if transitions.len() >= MAX_DECOY_TRANSITIONS {
                        return Err(Error::InvalidValue(format!(
                            "decoy generation exceeds {MAX_DECOY_TRANSITIONS} transitions"
                        )));
                    }
                    transitions.push(decoy_tr);
                } else {
                    exclusion.insert(decoy_tr.peptide_ref().to_owned());
                }
            }
        }
        transitions.retain(|tr| !exclusion.contains(tr.peptide_ref()));

        let mut protein_ids: HashSet<String> = HashSet::new();
        let decoy_peptides: Vec<Peptide> = peptides
            .into_iter()
            .filter(|p| !exclusion.contains(&p.id))
            .inspect(|p| protein_ids.extend(p.protein_refs.iter().cloned()))
            .collect();
        let decoy_proteins: Vec<Protein> = proteins
            .into_iter()
            .filter(|p| protein_ids.contains(&p.id))
            .collect();

        dec.set_transitions(transitions);
        dec.set_peptides(decoy_peptides);
        dec.set_proteins(decoy_proteins);
        Ok(())
    }
}

/// `getModifiedPeptideSequence_(peptide) + toStr(getChargeState())`; an unset
/// charge prints as `0`, the source release build's value.
fn precursor_key(peptide: &Peptide) -> String {
    format!(
        "{}{}",
        MRMDecoy::modified_peptide_sequence(peptide),
        peptide.charge_state().unwrap_or(0)
    )
}

/// Erase the fixed positions from `index`, back to front, as the source's
/// `vector::erase` loop (positions shift after each erasure).
fn remove_fixed(index: &mut Vec<usize>, fixed: &[usize]) -> Result<()> {
    for &position in fixed.iter().rev() {
        if position >= index.len() {
            return Err(fixed_range_error());
        }
        index.remove(position);
    }
    Ok(())
}

/// Re-insert each fixed position at its own index, front to back.
fn reinsert_fixed(index: &mut Vec<usize>, fixed: &[usize]) -> Result<()> {
    for &position in fixed {
        if position > index.len() {
            return Err(fixed_range_error());
        }
        index.insert(position, position);
    }
    Ok(())
}

fn fixed_range_error() -> Error {
    Error::InvalidValue(
        "fixed-residue positions run out of range (repeated letter in the constant pattern)".into(),
    )
}

/// Move every residue modification (`0 <= location < len`) to the new
/// position of its residue.
fn relocate(locations: &mut [i32], peptide_index: &[usize], len: usize) {
    for location in locations {
        let Ok(old) = usize::try_from(*location) else {
            continue;
        };
        if old >= len {
            continue;
        }
        if let Some(new) = peptide_index.iter().position(|&from| from == old) {
            *location = i32::try_from(new).unwrap_or(*location);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boost_uniform_int_ranges() {
        let mut rng = BoostUniformInt::new(42);
        for _ in 0..1000 {
            assert!(rng.draw() <= 9);
            assert!(rng.below(3) < 3);
        }
        assert_eq!(rng.below(1), 0);
        assert_eq!(rng.below(0), 0);
    }

    #[test]
    fn fixed_position_bookkeeping_rejects_out_of_range() {
        let mut index = vec![0, 1];
        assert!(remove_fixed(&mut index, &[5]).is_err());
        let mut index = vec![0];
        assert!(reinsert_fixed(&mut index, &[3]).is_err());
    }
}
