// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Fragment index for database search, ported from
//! `ANALYSIS/ID/FragmentIndex.h`.
//!
//! From a set of FASTA entries, the index stores the theoretical m/z of every
//! selected fragment ion (b and y by default) of every peptide generated from
//! the database, arranged in two dimensions: fragments are sorted by their own
//! m/z and cut into buckets, and within a bucket they are sorted by the
//! peptide they come from, peptides themselves being sorted by precursor
//! (M+H)+ mass. A query walks the buckets within the fragment tolerance of a
//! peak and binary-searches the precursor-compatible peptide range inside
//! each.
//!
//! The SNES ("Speedy Non-specific Enzyme Search") mode indexes mother
//! peptides for non-specific digestion; its build and realization are ported,
//! its spectrum query is not (see [`FragmentIndex::query_spectrum`]).
//!
//! Parameters: [`FragmentIndexParameters`]. See
//! `docs/FRAGMENT_INDEX_SUPPORT.md` for the API mapping and evidence.

#[path = "fragment_index_parameters.rs"]
mod parameters;

pub use parameters::FragmentIndexParameters;

use crate::chemistry::{
    AASequence, C13C12_MASSDIFF_U, EmpiricalFormula, ModifiedPeptideGenerator, PROTON_MASS_U,
    Protease, ProteaseDigestion, ResidueModification, TermSpecificity,
};
use crate::concept::parallel::{Threads, map_collect};
use crate::format::fasta::FASTAEntry;
use crate::{Error, MSSpectrum, Peak1D, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock};

/// Bit 31 of [`Peptide::mod_bitmask`] in SNES mode: set for a Single-C
/// (C-terminus anchored) mother, clear for a Single-N mother.
pub const SNES_KIND_BIT_MASK: u32 = 1 << 31;
/// Bits 0..30 of [`Peptide::mod_bitmask`] in SNES mode.
pub const SNES_SLOT_MASK: u32 = !SNES_KIND_BIT_MASK;
/// Maximum variable modification slots per peptide (a `u32` bitmask).
pub const MAX_MOD_SLOTS: usize = 32;
/// Longest FASTA entry the index accepts: peptide coordinates are 16-bit.
pub const MAX_PROTEIN_LENGTH: usize = u16::MAX as usize;
/// Fragments per bucket, chosen by the source to approximate a 0.02 Da
/// fragment-bin density.
pub const BUCKET_SIZE: usize = 4096;

/// Compact descriptor of a peptide held by the index (source `Peptide`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Peptide {
    /// Index into the FASTA entries passed to [`FragmentIndex::build`].
    pub protein_idx: u32,
    /// Active variable modification slots (bit `s` = slot `s` of
    /// [`FragmentIndex::build_mod_slots`]); 0 = unmodified or fixed-only, up
    /// to 32 slots. In SNES mode bit 31 is [`SNES_KIND_BIT_MASK`].
    pub mod_bitmask: u32,
    /// `(start, length)` of the peptide within the protein, 0-based start,
    /// length in residues (as `std::string::substr(start, length)`).
    pub sequence: (u16, u16),
    /// Monoisotopic (M+H)+ of this peptide; the sort key of the index.
    pub precursor_mz: f32,
}

impl Peptide {
    /// Construct a descriptor (source constructor).
    pub fn new(
        protein_idx: u32,
        mod_bitmask: u32,
        sequence: (u16, u16),
        precursor_mz: f32,
    ) -> Self {
        Self {
            protein_idx,
            mod_bitmask,
            sequence,
            precursor_mz,
        }
    }
}

/// A match between a query spectrum and an indexed peptide (source
/// `SpectrumMatch`).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpectrumMatch {
    /// Number of peak-fragment hits.
    pub num_matched: u32,
    /// SNES: active slots of the realized sub-peptide; 0 = unmodified. Unused
    /// outside SNES mode.
    pub subset_bitmask: u32,
    /// SNES: sum of the variable modification deltas of this match; 0
    /// outside SNES mode.
    pub sigma_delta: f32,
    /// The precursor charge the search used.
    pub precursor_charge: u16,
    /// The isotope error the search used.
    pub isotope_error: i16,
    /// Index of the matched peptide in [`FragmentIndex::peptides`].
    pub peptide_idx: usize,
}

/// Container of [`SpectrumMatch`]es (source `SpectrumMatchesTopN`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SpectrumMatchesTopN {
    /// The preliminary candidates.
    pub hits: Vec<SpectrumMatch>,
}

impl SpectrumMatchesTopN {
    /// Remove every candidate.
    pub fn clear(&mut self) {
        self.hits.clear();
    }
}

impl std::ops::AddAssign<&SpectrumMatchesTopN> for SpectrumMatchesTopN {
    /// Append the other container's candidates (source `operator+=`).
    fn add_assign(&mut self, other: &SpectrumMatchesTopN) {
        self.hits.extend_from_slice(&other.hits);
    }
}

/// A match between one query peak and one indexed fragment (source `Hit`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    /// Index of the peptide in [`FragmentIndex::peptides`].
    pub peptide_idx: u32,
    /// The fragment's m/z.
    pub fragment_mz: f32,
}

/// One entry of the index (source protected `Fragment`), exposed read-only
/// through [`FragmentIndex::fragments`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Fragment {
    /// Index of the originating peptide.
    pub peptide_idx: u32,
    /// Fragment m/z at charge 1.
    pub fragment_mz: f32,
}

/// Protein-anchor restriction of an SNES bin walk (source `SnesAnchor`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SnesAnchor {
    /// No anchor restriction.
    #[default]
    None,
    /// The mother must start at protein position 0.
    ProtNterm,
    /// The mother must end at the protein's last residue.
    ProtCterm,
}

/// A variable modification slot of one peptide (source protected `ModSlot`).
#[derive(Clone, Debug, PartialEq)]
pub struct ModSlot {
    /// Residue index, or [`ModSlot::NTERM_SLOT`] / [`ModSlot::CTERM_SLOT`].
    pub position: u16,
    /// The modification's monoisotopic mass delta.
    pub delta_mass: f64,
    /// The modification itself, for sequence reconstruction.
    pub modification: Arc<ResidueModification>,
}

impl ModSlot {
    /// Position sentinel of a pure N-terminal modification.
    pub const NTERM_SLOT: u16 = u16::MAX - 1;
    /// Position sentinel of a pure C-terminal modification.
    pub const CTERM_SLOT: u16 = u16::MAX;
}

/// A variable modification in the per-residue lookup (source `VarModEntry`).
#[derive(Clone, Debug)]
struct VarModEntry {
    delta_mass: f64,
    modification: Arc<ResidueModification>,
    term_spec: TermSpecificity,
}

/// Modification lookup tables (source `initModificationTables_` members).
#[derive(Clone, Debug)]
struct ModTables {
    fixed_deltas: [f64; 128],
    fixed_mods: Vec<Option<Arc<ResidueModification>>>,
    fixed_nterm_delta: f64,
    fixed_cterm_delta: f64,
    fixed_nterm_mod: Option<Arc<ResidueModification>>,
    fixed_cterm_mod: Option<Arc<ResidueModification>>,
    variable: Vec<Vec<VarModEntry>>,
    variable_nterm: Vec<VarModEntry>,
    variable_cterm: Vec<VarModEntry>,
}

/// Residue internal masses and ion offsets (source static tables).
struct MassTables {
    residue: [f64; 128],
    water: f64,
    b: f64,
    y: f64,
    a: f64,
    c: f64,
    x: f64,
    z: f64,
}

fn mass_tables() -> Result<&'static MassTables> {
    static TABLES: OnceLock<std::result::Result<MassTables, String>> = OnceLock::new();
    TABLES
        .get_or_init(|| build_mass_tables().map_err(|e| e.to_string()))
        .as_ref()
        .map_err(|e| Error::InvalidValue(format!("residue mass table: {e}")))
}

fn build_mass_tables() -> Result<MassTables> {
    let formula = |text: &str| -> Result<f64> { Ok(EmpiricalFormula::parse(text)?.mono_mass()) };
    let water = formula("H2O")?;
    let mut residue = [0.0; 128];
    for letter in b'A'..=b'Z' {
        // Residues without a known composition (B, Z, X) keep 0, which the
        // index never reads: sequences containing them are skipped.
        let text = char::from(letter).to_string();
        if let Ok(peptide) = AASequence::parse(&text) {
            if let Ok(full) = peptide.formula() {
                residue[usize::from(letter)] = full.mono_mass() - water;
            }
        }
    }
    // Residue::getInternalTo*Ion() at the pinned source: b = H - H, y = OH + H,
    // a = H - CHO, c = H + NH2, x = OH + CO - H, z = OH - NH2.
    Ok(MassTables {
        residue,
        water,
        b: 0.0,
        y: formula("H2O")?,
        a: formula("C-1O-1")?,
        c: formula("H3N")?,
        x: formula("CO2")?,
        z: formula("H-1N-1O")?,
    })
}

/// Which ion series a fragment generation emits.
#[derive(Clone, Copy, Debug)]
struct Series {
    b: bool,
    a: bool,
    c: bool,
    y: bool,
    x: bool,
    z: bool,
}

/// Fragment index (source `FragmentIndex`).
///
/// Construct with [`FragmentIndex::new`] or [`FragmentIndex::from_param`],
/// then [`FragmentIndex::build`] before querying. The source's
/// `DefaultParamHandler` base is [`FragmentIndexParameters`].
#[derive(Clone, Debug)]
pub struct FragmentIndex {
    params: FragmentIndexParameters,
    tables: ModTables,
    sigma_delta_set: Vec<f64>,
    sigma_delta_set_with_prot_nterm: Vec<f64>,
    sigma_delta_set_with_prot_cterm: Vec<f64>,
    peptides: Vec<Peptide>,
    fragments: Vec<Fragment>,
    protein_lengths: Vec<u32>,
    bucket_min_mz: Vec<f32>,
    bucket_size: usize,
    is_build: bool,
    /// Worker threads for building (native; the source uses OpenMP's
    /// setting). The index is identical for every value.
    pub threads: Threads,
}

impl FragmentIndex {
    /// An unbuilt index with these parameters (source constructor followed by
    /// `setParameters`, whose `updateMembers_` builds the modification tables
    /// and the SNES Σ sets).
    ///
    /// Logs a warning when the precursor window triggers open-search mode, and
    /// when an SNES Σ set exceeds 64 entries.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for an invalid parameter combination
    /// (see [`FragmentIndexParameters::validate`]) and for an unknown
    /// modification name.
    pub fn new(params: FragmentIndexParameters) -> Result<Self> {
        params.validate()?;
        mass_tables()?;
        if params.is_open_search_mode() {
            let _ = crate::openms_log_warn!(
                "[FragmentIndex] Open-search mode auto-triggered: window [-{}, +{}] {} exceeds threshold. Isotope-error iteration collapses to [0, 0].",
                params.precursor_mass_tolerance_lower,
                params.precursor_mass_tolerance_upper,
                if params.precursor_mass_tolerance_unit_ppm {
                    "ppm"
                } else {
                    "Da"
                }
            );
        }
        let tables = init_modification_tables(&params)?;
        let mut index = Self {
            params,
            tables,
            sigma_delta_set: Vec::new(),
            sigma_delta_set_with_prot_nterm: Vec::new(),
            sigma_delta_set_with_prot_cterm: Vec::new(),
            peptides: Vec::new(),
            fragments: Vec::new(),
            protein_lengths: Vec::new(),
            bucket_min_mz: Vec::new(),
            bucket_size: BUCKET_SIZE,
            is_build: false,
            threads: Threads::default(),
        };
        index.sigma_delta_set = index.compute_snes_sigma_delta_set(false, false);
        index.sigma_delta_set_with_prot_nterm = index.compute_snes_sigma_delta_set(true, false);
        index.sigma_delta_set_with_prot_cterm = index.compute_snes_sigma_delta_set(false, true);
        let largest = index
            .sigma_delta_set
            .len()
            .max(index.sigma_delta_set_with_prot_nterm.len())
            .max(index.sigma_delta_set_with_prot_cterm.len());
        if index.is_snes_mode() && largest > 64 {
            let _ = crate::openms_log_warn!(
                "[FragmentIndex] SNES Σ_delta set has {largest} entries — query performance will scale linearly with this. Consider reducing modifications:variable or variable_max_per_peptide."
            );
        }
        Ok(index)
    }

    /// An unbuilt index from a [`crate::param::Param`], as the source's
    /// `setParameters`.
    ///
    /// # Errors
    ///
    /// As [`FragmentIndexParameters::from_param`] and [`FragmentIndex::new`].
    pub fn from_param(parameters: &crate::param::Param) -> Result<Self> {
        Self::new(FragmentIndexParameters::from_param(parameters)?)
    }

    /// Replace the parameters (source `setParameters`).
    ///
    /// The source keeps an existing index when parameters change, although
    /// its peptides' modification bitmasks then refer to the old
    /// modification tables; this port clears the index instead, so a rebuild
    /// is required. On error the index is unchanged.
    ///
    /// # Errors
    ///
    /// As [`FragmentIndex::new`].
    pub fn set_parameters(&mut self, params: FragmentIndexParameters) -> Result<()> {
        let threads = self.threads;
        *self = Self::new(params)?;
        self.threads = threads;
        Ok(())
    }

    /// The current parameters.
    pub fn parameters(&self) -> &FragmentIndexParameters {
        &self.params
    }

    /// Whether [`FragmentIndex::build`] has completed (source `isBuild`).
    pub fn is_build(&self) -> bool {
        self.is_build
    }

    /// All peptides of the index, sorted by precursor m/z (source
    /// `getPeptides`); empty before [`FragmentIndex::build`].
    pub fn peptides(&self) -> &[Peptide] {
        &self.peptides
    }

    /// Number of fragments stored in the index (source `getNumFragments`); 0
    /// before [`FragmentIndex::build`].
    pub fn num_fragments(&self) -> usize {
        self.fragments.len()
    }

    /// The fragment array: sorted by m/z into buckets of
    /// [`FragmentIndex::bucket_size`], each bucket sorted by peptide index
    /// (the source's protected `fi_fragments_`, exposed read-only).
    pub fn fragments(&self) -> &[Fragment] {
        &self.fragments
    }

    /// Fragments per bucket: [`BUCKET_SIZE`], or 1 for an empty index.
    pub fn bucket_size(&self) -> usize {
        self.bucket_size
    }

    /// Whether the index was built in SNES mode: `snes_enabled` together with
    /// `peptide:enzyme_specificity` `none` (source `isSnesMode`).
    pub fn is_snes_mode(&self) -> bool {
        self.params.is_snes_mode()
    }

    /// Open-search auto-detection: `max(lower, upper)` strictly above 1000
    /// ppm, or 1 Da (source static `isOpenSearchMode`, shared with ProSE and
    /// the TOPP tool).
    pub fn is_open_search_mode(lower_magnitude: f64, upper_magnitude: f64, unit_ppm: bool) -> bool {
        let threshold = if unit_ppm { 1000.0 } else { 1.0 };
        lower_magnitude.max(upper_magnitude) > threshold
    }

    /// Whether a bitmask describes a Single-C mother (SNES only).
    pub fn is_single_c_mother(mod_bitmask: u32) -> bool {
        mod_bitmask & SNES_KIND_BIT_MASK != 0
    }

    /// Whether a bitmask describes a Single-N mother (SNES only).
    pub fn is_single_n_mother(mod_bitmask: u32) -> bool {
        mod_bitmask & SNES_KIND_BIT_MASK == 0
    }

    /// The SNES Σ-delta set for a bin-walk anchor: baseline, with
    /// protein-N-terminal-only, or with protein-C-terminal-only variable
    /// modifications (the source's protected `snes_sigma_delta_set_*`).
    pub fn snes_sigma_delta_set(&self, anchor: SnesAnchor) -> &[f64] {
        match anchor {
            SnesAnchor::None => &self.sigma_delta_set,
            SnesAnchor::ProtNterm => &self.sigma_delta_set_with_prot_nterm,
            SnesAnchor::ProtCterm => &self.sigma_delta_set_with_prot_cterm,
        }
    }

    /// Delete the index (source `clear`); [`FragmentIndex::is_build`] becomes
    /// `false`. Parameters and modification tables are kept.
    pub fn clear(&mut self) {
        self.fragments = Vec::new();
        self.peptides = Vec::new();
        self.bucket_min_mz = Vec::new();
        self.protein_lengths = Vec::new();
        self.is_build = false;
    }

    /// Distinct sums of variable modification deltas reachable with at most
    /// `modifications:variable_max_per_peptide` modifications, ascending and
    /// always including 0 (source protected `computeSnesSigmaDeltaSet_`).
    ///
    /// Configuration-global: residue inventory is not considered. Deltas
    /// within 1e-6 Da are one value, both among the configured deltas and
    /// among the sums. Protein-terminal-only modifications participate only
    /// when requested.
    pub fn compute_snes_sigma_delta_set(
        &self,
        include_prot_nterm_mods: bool,
        include_prot_cterm_mods: bool,
    ) -> Vec<f64> {
        let tables = &self.tables;
        let mut eligible = Vec::new();
        let mut collect = |entries: &[VarModEntry]| {
            for e in entries {
                if e.term_spec == TermSpecificity::ProteinNTerm && !include_prot_nterm_mods {
                    continue;
                }
                if e.term_spec == TermSpecificity::ProteinCTerm && !include_prot_cterm_mods {
                    continue;
                }
                eligible.push(e.delta_mass);
            }
        };
        collect(&tables.variable_nterm);
        collect(&tables.variable_cterm);
        for per_aa in &tables.variable {
            collect(per_aa);
        }
        eligible.sort_by(f64::total_cmp);
        eligible.dedup_by(|b, a| (*a - *b).abs() < 1e-6);
        let mut result = vec![0.0];
        let max = self.params.modifications_max_variable_mods_per_peptide;
        if eligible.is_empty() || max == 0 {
            return result;
        }
        let mut previous = vec![0.0];
        for _ in 0..max {
            let mut next = Vec::with_capacity(previous.len() * eligible.len());
            for &p in &previous {
                for &d in &eligible {
                    next.push(p + d);
                }
            }
            next.sort_by(f64::total_cmp);
            next.dedup_by(|b, a| (*a - *b).abs() < 1e-6);
            for &v in &next {
                let at = result.partition_point(|&r| r < v - 1e-6);
                if at == result.len() || (result[at] - v).abs() >= 1e-6 {
                    result.insert(at, v);
                }
            }
            previous = next;
        }
        result
    }

    /// The variable modification slots of a peptide (source protected
    /// `buildModSlots_`), at most [`MAX_MOD_SLOTS`].
    ///
    /// Deterministic order: pure N-terminal modifications, then residue
    /// modifications left to right (anywhere, plus peptide- or
    /// protein-terminal ones at the matching end), then pure C-terminal
    /// modifications. Protein-terminal modifications need
    /// `is_protein_nterm`/`is_protein_cterm`. Slots beyond 32 are dropped
    /// silently, as in the source.
    pub fn build_mod_slots(
        &self,
        sequence: &[u8],
        is_protein_nterm: bool,
        is_protein_cterm: bool,
    ) -> Vec<ModSlot> {
        let tables = &self.tables;
        let mut slots = Vec::new();
        for e in &tables.variable_nterm {
            if e.term_spec == TermSpecificity::ProteinNTerm && !is_protein_nterm {
                continue;
            }
            if slots.len() >= MAX_MOD_SLOTS {
                break;
            }
            slots.push(slot(ModSlot::NTERM_SLOT, e));
        }
        let last = sequence.len().wrapping_sub(1);
        for (i, &aa) in sequence.iter().enumerate() {
            for e in tables
                .variable
                .get(usize::from(aa))
                .map_or(&[][..], Vec::as_slice)
            {
                if slots.len() >= MAX_MOD_SLOTS {
                    break;
                }
                let applies = match e.term_spec {
                    TermSpecificity::Anywhere => true,
                    TermSpecificity::NTerm => i == 0,
                    TermSpecificity::ProteinNTerm => i == 0 && is_protein_nterm,
                    TermSpecificity::CTerm => i == last,
                    TermSpecificity::ProteinCTerm => i == last && is_protein_cterm,
                };
                if applies {
                    // Positions fit u16: proteins are at most 65535 residues.
                    slots.push(slot(i as u16, e));
                }
            }
        }
        for e in &tables.variable_cterm {
            if e.term_spec == TermSpecificity::ProteinCTerm && !is_protein_cterm {
                continue;
            }
            if slots.len() >= MAX_MOD_SLOTS {
                break;
            }
            slots.push(slot(ModSlot::CTERM_SLOT, e));
        }
        slots
    }

    /// Build the index from FASTA entries (source `build`).
    ///
    /// Generates peptides (enzymatic digestion, or SNES mothers), their
    /// fragments within `fragment:min_mz`..`fragment:max_mz` (skipping ions
    /// with index at most `fragment:min_ion_index`), sorts fragments by m/z
    /// (then peptide index), cuts them into buckets of [`BUCKET_SIZE`]
    /// recording each bucket's smallest m/z, and sorts each bucket by peptide
    /// index (stable, as the source's radix sort). Peptides are sorted by
    /// `(precursor_mz, protein_idx)` with a **stable** sort, so peptides with
    /// equal keys keep generation order; the source's `std::sort` leaves that
    /// order to the library. An empty result is a built, empty index.
    ///
    /// A previous index is replaced.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a FASTA entry longer than
    /// [`MAX_PROTEIN_LENGTH`] (the source's `Exception::InvalidParameter`),
    /// or with characters other than `A`-`Z` (the source indexes a 128-entry
    /// table with them, undefined above 127), and propagates digestion errors
    /// and resource ceilings (digestion products and work per protein,
    /// [`MAX_MOD_SLOTS`] enumeration). On error the index is cleared.
    pub fn build(&mut self, fasta_entries: &[FASTAEntry]) -> Result<()> {
        self.clear();
        let mut protein_lengths = Vec::with_capacity(fasta_entries.len());
        for entry in fasta_entries {
            if entry.sequence.len() > MAX_PROTEIN_LENGTH {
                return Err(Error::InvalidValue(format!(
                    "FragmentIndex: FASTA entry '{}' has {} residues, exceeding the supported maximum of 65535 (peptide offsets are stored as 16-bit). Split long contigs / six-frame-translated frames into windows of at most 65535 residues (with overlap >= peptide:max_size so no peptide is lost across a split) before building the index.",
                    entry.identifier,
                    entry.sequence.len()
                )));
            }
            if !entry.sequence.bytes().all(|b| b.is_ascii_uppercase()) {
                return Err(Error::InvalidValue(format!(
                    "FragmentIndex: FASTA entry '{}' contains characters other than A-Z",
                    entry.identifier
                )));
            }
            protein_lengths.push(entry.sequence.len() as u32);
        }
        let peptides = if self.is_snes_mode() {
            self.generate_snes_mothers(fasta_entries)?
        } else {
            self.generate_peptides(fasta_entries)?
        };
        let chunks: Vec<std::ops::Range<usize>> = (0..peptides.len())
            .step_by(4096)
            .map(|s| s..(s + 4096).min(peptides.len()))
            .collect();
        let generated = map_collect(&chunks, self.threads, |range| {
            let mut out = Vec::new();
            for idx in range.clone() {
                self.fragments_for_peptide(&mut out, &peptides[idx], idx as u32, fasta_entries)?;
            }
            Ok::<_, Error>(out)
        });
        let mut fragments = Vec::new();
        for chunk in generated {
            fragments.extend(chunk?);
        }
        fragments.sort_by(|a: &Fragment, b: &Fragment| {
            a.fragment_mz
                .total_cmp(&b.fragment_mz)
                .then(a.peptide_idx.cmp(&b.peptide_idx))
        });
        self.peptides = peptides;
        self.protein_lengths = protein_lengths;
        if fragments.is_empty() {
            self.bucket_size = 1;
            let _ =
                crate::openms_log_info!("[FragmentIndex] No fragments generated — index is empty.");
            self.is_build = true;
            return Ok(());
        }
        self.bucket_size = BUCKET_SIZE;
        let mut bucket_min_mz = Vec::with_capacity(fragments.len().div_ceil(BUCKET_SIZE));
        for bucket in fragments.chunks_mut(BUCKET_SIZE) {
            bucket_min_mz.push(bucket[0].fragment_mz);
            bucket.sort_by_key(|f| f.peptide_idx);
        }
        self.fragments = fragments;
        self.bucket_min_mz = bucket_min_mz;
        self.is_build = true;
        Ok(())
    }

    /// Enzymatic peptides with their fixed/variable modification variants
    /// (source protected `generatePeptides`, non-SNES branch).
    fn generate_peptides(&self, fasta_entries: &[FASTAEntry]) -> Result<Vec<Peptide>> {
        let p = &self.params;
        let masses = mass_tables()?;
        let has_variable_mods = !p.modifications_variable.is_empty();
        let digestor = ProteaseDigestion {
            enzyme: Protease::from_name(&p.enzyme)?,
            missed_cleavages: p.missed_cleavages,
            min_length: p.peptide_min_length,
            max_length: (p.peptide_max_length > 0).then_some(p.peptide_max_length),
            specificity: p.enzyme_specificity,
            ..ProteaseDigestion::default()
        };
        let tables = &self.tables;
        let per_protein = map_collect(fasta_entries, self.threads, |protein| {
            let mut out = Vec::new();
            let mut skipped = 0_usize;
            let sequence = protein.sequence.as_bytes();
            for product in digestor.digest_ranges(&protein.sequence)? {
                let (start, len) = (product.start, product.end - product.start);
                let sub = &sequence[start..start + len];
                if sub.iter().any(|b| matches!(b, b'X' | b'B' | b'Z')) {
                    skipped += 1;
                    continue;
                }
                let mut base_mass = masses.water
                    + PROTON_MASS_U
                    + tables.fixed_nterm_delta
                    + tables.fixed_cterm_delta;
                for &aa in sub {
                    base_mass +=
                        masses.residue[usize::from(aa)] + tables.fixed_deltas[usize::from(aa)];
                }
                let coordinates = (start as u16, len as u16);
                let mut push = |bitmask: u32, mass: f64| {
                    let mz = mass as f32;
                    if p.peptide_min_mass as f32 <= mz && mz <= p.peptide_max_mass as f32 {
                        out.push((bitmask, coordinates, mz));
                    }
                };
                if has_variable_mods {
                    let slots =
                        self.build_mod_slots(sub, start == 0, start + len == sequence.len());
                    if slots.is_empty() {
                        push(0, base_mass);
                    } else {
                        let n = slots.len();
                        let mut conflict = [0_u32; MAX_MOD_SLOTS];
                        for a in 0..n {
                            for b in a + 1..n {
                                if slots[a].position == slots[b].position {
                                    conflict[a] |= 1 << b;
                                    conflict[b] |= 1 << a;
                                }
                            }
                        }
                        // n <= 32; the source's `1u << 32` is undefined for 32 slots.
                        let max_bitmask: u64 = 1_u64 << n;
                        for bitmask in 0..max_bitmask {
                            let bitmask = bitmask as u32;
                            if bitmask.count_ones() as usize
                                > p.modifications_max_variable_mods_per_peptide
                            {
                                continue;
                            }
                            let conflicting = (0..n).any(|s| {
                                bitmask & (1 << s) != 0 && bitmask & conflict[s] & !(1 << s) != 0
                            });
                            if conflicting {
                                continue;
                            }
                            let mut variant = base_mass;
                            for (s, slot) in slots.iter().enumerate() {
                                if bitmask & (1 << s) != 0 {
                                    variant += slot.delta_mass;
                                }
                            }
                            push(bitmask, variant);
                        }
                    }
                } else {
                    push(0, base_mass);
                }
            }
            Ok::<_, Error>((out, skipped))
        });
        let mut peptides = Vec::new();
        let mut skipped = 0;
        for (protein_idx, result) in per_protein.into_iter().enumerate() {
            let (list, s) = result?;
            skipped += s;
            peptides.extend(list.into_iter().map(|(bitmask, coordinates, mz)| {
                Peptide::new(protein_idx as u32, bitmask, coordinates, mz)
            }));
        }
        if skipped > 0 {
            let _ = crate::openms_log_warn!(
                "{skipped} peptides skipped due to unknown or ambiguous AA (X/B/Z)"
            );
        }
        sort_peptides(&mut peptides);
        Ok(peptides)
    }

    /// SNES mother peptides (source protected `generateSNESMothers_`).
    ///
    /// For each maximal span of a protein without `X`/`B`/`Z`: a Single-N
    /// mother at every anchor `i` with at least `peptide:min_size` residues
    /// left, spanning `min(max_size, span_end - i)` residues (`max_size` 0 =
    /// no limit), and a Single-C mother ending at every `j`, spanning back up
    /// to `max_size` residues. Mothers lighter than `peptide:min_mass` are
    /// dropped; the upper mass bound is left to realization, as in the
    /// source. Variable modifications are not enumerated on mothers.
    fn generate_snes_mothers(&self, fasta_entries: &[FASTAEntry]) -> Result<Vec<Peptide>> {
        let p = &self.params;
        let masses = mass_tables()?;
        let tables = &self.tables;
        if !p.modifications_variable.is_empty() {
            let _ = crate::openms_log_info!(
                "[FragmentIndex] SNES v1.1: variable modifications are not enumerated on mother peptides (mother index stores unmodified mothers); query-time subset enumeration applies them to realized sub-peptides using configured mods: {}",
                p.modifications_variable.join(", ")
            );
        }
        let base_constants =
            masses.water + PROTON_MASS_U + tables.fixed_nterm_delta + tables.fixed_cterm_delta;
        let min_length = p.peptide_min_length;
        let snes_min_length = min_length.max(1);
        let per_protein = map_collect(fasta_entries, self.threads, |protein| {
            let seq = protein.sequence.as_bytes();
            let l = seq.len();
            let mut out = Vec::new();
            let mut skipped = 0_usize;
            if l < min_length {
                return (out, skipped);
            }
            let max_length = if p.peptide_max_length == 0 {
                l
            } else {
                p.peptide_max_length
            };
            let mut emit = |start: usize, length: usize, single_c: bool| {
                if length < min_length {
                    return;
                }
                let mut mass = base_constants;
                for &aa in &seq[start..start + length] {
                    mass += masses.residue[usize::from(aa)] + tables.fixed_deltas[usize::from(aa)];
                }
                let mz = mass as f32;
                if mz < p.peptide_min_mass as f32 {
                    return;
                }
                let kind = if single_c { SNES_KIND_BIT_MASK } else { 0 };
                out.push((kind, (start as u16, length as u16), mz));
            };
            let mut sweep = |s: usize, e: usize, emit: &mut dyn FnMut(usize, usize, bool)| {
                if e <= s || e - s < min_length {
                    if e > s {
                        skipped += 1;
                    }
                    return;
                }
                let mut i = s;
                while i + min_length <= e {
                    emit(i, max_length.min(e - i), false);
                    i += 1;
                }
                for j in (s + snes_min_length - 1)..e {
                    let length = max_length.min(j + 1 - s);
                    emit(j + 1 - length, length, true);
                }
            };
            let mut start = 0;
            let mut bad = seq.iter().position(|b| matches!(b, b'X' | b'B' | b'Z'));
            loop {
                match bad {
                    None => {
                        sweep(start, l, &mut emit);
                        break;
                    }
                    Some(at) => {
                        sweep(start, at, &mut emit);
                        start = at + 1;
                        if start >= l {
                            break;
                        }
                        bad = seq[start..]
                            .iter()
                            .position(|b| matches!(b, b'X' | b'B' | b'Z'))
                            .map(|offset| start + offset);
                    }
                }
            }
            (out, skipped)
        });
        let mut peptides = Vec::new();
        let mut skipped = 0;
        for (protein_idx, (list, s)) in per_protein.into_iter().enumerate() {
            skipped += s;
            peptides.extend(list.into_iter().map(|(bitmask, coordinates, mz)| {
                Peptide::new(protein_idx as u32, bitmask, coordinates, mz)
            }));
        }
        sort_peptides(&mut peptides);
        let _ = crate::openms_log_info!(
            "Generated {} SNES mothers ({skipped} spans skipped — shorter than peptide:min_size).",
            peptides.len()
        );
        Ok(peptides)
    }

    /// The fragments of one peptide (the per-peptide body of source `build`).
    fn fragments_for_peptide(
        &self,
        out: &mut Vec<Fragment>,
        pep: &Peptide,
        peptide_idx: u32,
        fasta_entries: &[FASTAEntry],
    ) -> Result<()> {
        let p = &self.params;
        let tables = &self.tables;
        let protein = fasta_entries
            .get(pep.protein_idx as usize)
            .ok_or_else(|| Error::InvalidValue("peptide protein index out of range".into()))?;
        let (start, len) = (usize::from(pep.sequence.0), usize::from(pep.sequence.1));
        let seq = protein
            .sequence
            .as_bytes()
            .get(start..start + len)
            .ok_or_else(|| Error::InvalidValue("peptide coordinates out of range".into()))?;
        let fixed = |seq: &[u8]| -> (Vec<f64>, bool) {
            let mut masses = vec![0.0; seq.len()];
            let mut any = false;
            for (i, &aa) in seq.iter().enumerate() {
                let delta = tables.fixed_deltas[usize::from(aa)];
                if delta != 0.0 {
                    masses[i] = delta;
                    any = true;
                }
            }
            (masses, any)
        };
        let has_modifications =
            !(p.modifications_fixed.is_empty() && p.modifications_variable.is_empty());
        if self.is_snes_mode() {
            let single_c = Self::is_single_c_mother(pep.mod_bitmask);
            let (masses, any) = fixed(seq);
            let n_term = if single_c {
                0.0
            } else {
                tables.fixed_nterm_delta
            };
            let c_term = if single_c {
                tables.fixed_cterm_delta
            } else {
                0.0
            };
            let series = Series {
                b: !single_c,
                a: false,
                c: false,
                y: single_c,
                x: false,
                z: false,
            };
            self.generate_fragments_for_series(
                out,
                seq,
                peptide_idx,
                n_term,
                c_term,
                any.then_some(&masses[..]),
                series,
            )?;
        } else if !has_modifications {
            self.generate_fragments_for_series(
                out,
                seq,
                peptide_idx,
                0.0,
                0.0,
                None,
                self.class_series(),
            )?;
        } else if pep.mod_bitmask == 0 {
            let (masses, any) = fixed(seq);
            self.generate_fragments_for_series(
                out,
                seq,
                peptide_idx,
                tables.fixed_nterm_delta,
                tables.fixed_cterm_delta,
                any.then_some(&masses[..]),
                self.class_series(),
            )?;
        } else {
            let (mut masses, mut any) = fixed(seq);
            let (mut n_term, mut c_term) = (tables.fixed_nterm_delta, tables.fixed_cterm_delta);
            let slot_bits = pep.mod_bitmask & SNES_SLOT_MASK;
            let slots =
                self.build_mod_slots(seq, start == 0, start + len == protein.sequence.len());
            for (s, slot) in slots.iter().enumerate() {
                if slot_bits & (1 << s) == 0 {
                    continue;
                }
                match slot.position {
                    ModSlot::NTERM_SLOT => n_term += slot.delta_mass,
                    ModSlot::CTERM_SLOT => c_term += slot.delta_mass,
                    position => {
                        masses[usize::from(position)] += slot.delta_mass;
                        any = true;
                    }
                }
            }
            self.generate_fragments_for_series(
                out,
                seq,
                peptide_idx,
                n_term,
                c_term,
                any.then_some(&masses[..]),
                self.class_series(),
            )?;
        }
        Ok(())
    }

    fn class_series(&self) -> Series {
        let p = &self.params;
        Series {
            b: p.add_b_ions,
            a: p.add_a_ions,
            c: p.add_c_ions,
            y: p.add_y_ions,
            x: p.add_x_ions,
            z: p.add_z_ions,
        }
    }

    /// Prefix and suffix fragment m/z at charge 1 from cumulative residue
    /// masses (source protected `generateFragmentsForSeries_`); masses are
    /// summed in `f64` and narrowed to `f32` per fragment.
    #[allow(clippy::too_many_arguments)]
    fn generate_fragments_for_series(
        &self,
        out: &mut Vec<Fragment>,
        seq: &[u8],
        peptide_idx: u32,
        n_term_mod_mass: f64,
        c_term_mod_mass: f64,
        residue_mod_masses: Option<&[f64]>,
        series: Series,
    ) -> Result<()> {
        let masses = mass_tables()?;
        let p = &self.params;
        let min_ion_index = p.min_ion_index;
        let mut emit = |cumulative: f64, offset: f64| {
            let mz = (cumulative + offset) as f32;
            if mz >= p.fragment_min_mz as f32 && mz <= p.fragment_max_mz as f32 {
                out.push(Fragment {
                    peptide_idx,
                    fragment_mz: mz,
                });
            }
        };
        if series.b || series.a || series.c {
            let mut cumulative = PROTON_MASS_U + n_term_mod_mass;
            for i in 0..seq.len().saturating_sub(1) {
                let mut residue = masses.residue[usize::from(seq[i])];
                if let Some(mods) = residue_mod_masses {
                    residue += mods[i];
                }
                cumulative += residue;
                if i + 1 <= min_ion_index {
                    continue;
                }
                if series.b {
                    emit(cumulative, masses.b);
                }
                if series.a {
                    emit(cumulative, masses.a);
                }
                if series.c {
                    emit(cumulative, masses.c);
                }
            }
        }
        if series.y || series.x || series.z {
            let mut cumulative = PROTON_MASS_U + c_term_mod_mass;
            let mut suffix_ion_num = 0;
            let mut j = seq.len();
            while j > 1 {
                let mut residue = masses.residue[usize::from(seq[j - 1])];
                if let Some(mods) = residue_mod_masses {
                    residue += mods[j - 1];
                }
                cumulative += residue;
                suffix_ion_num += 1;
                j -= 1;
                if suffix_ion_num <= min_ion_index {
                    continue;
                }
                if series.y {
                    emit(cumulative, masses.y);
                }
                if series.x {
                    emit(cumulative, masses.x);
                }
                if series.z {
                    emit(cumulative, masses.z);
                }
            }
        }
        Ok(())
    }

    /// The fully modified sequence of a peptide descriptor (source
    /// `reconstructModifiedSequence`, used for result output).
    ///
    /// Applies fixed residue modifications, fixed terminal modifications,
    /// then the variable modifications of the bitmask's slots (bit 31 masked
    /// off). A fixed terminal modification configured with a residue origin
    /// is applied to every peptide terminus whatever its residue, as the
    /// source does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a protein index or coordinates
    /// outside `fasta_entries`, and sequence errors.
    pub fn reconstruct_modified_sequence(
        &self,
        peptide: &Peptide,
        fasta_entries: &[FASTAEntry],
    ) -> Result<AASequence> {
        let (protein, start, len) =
            coordinates(peptide, fasta_entries, usize::from(peptide.sequence.1))?;
        let slot_bits = peptide.mod_bitmask & SNES_SLOT_MASK;
        self.sequence_with_mods(protein, start, len, slot_bits, true)
    }

    /// Shared body of the reconstruction functions: substring, fixed
    /// modifications (when any modification is configured), then the
    /// variable modifications of `slot_bits`.
    fn sequence_with_mods(
        &self,
        protein: &str,
        start: usize,
        len: usize,
        slot_bits: u32,
        require_modifications: bool,
    ) -> Result<AASequence> {
        let p = &self.params;
        let seq = AASequence::parse(&protein[start..start + len])?;
        let has_modifications =
            !(p.modifications_fixed.is_empty() && p.modifications_variable.is_empty());
        if !has_modifications && (require_modifications || slot_bits == 0) {
            return Ok(seq);
        }
        let tables = &self.tables;
        let bytes = &protein.as_bytes()[start..start + len];
        let mut residues: BTreeMap<usize, Arc<ResidueModification>> = BTreeMap::new();
        for (i, &aa) in bytes.iter().enumerate() {
            if let Some(Some(m)) = tables.fixed_mods.get(usize::from(aa)) {
                residues.insert(i, Arc::clone(m));
            }
        }
        let mut n_term = tables.fixed_nterm_mod.clone();
        let mut c_term = tables.fixed_cterm_mod.clone();
        if slot_bits != 0 {
            let slots = self.build_mod_slots(bytes, start == 0, start + len == protein.len());
            for (s, slot) in slots.iter().enumerate() {
                if slot_bits & (1 << s) == 0 {
                    continue;
                }
                match slot.position {
                    ModSlot::NTERM_SLOT => n_term = Some(Arc::clone(&slot.modification)),
                    ModSlot::CTERM_SLOT => c_term = Some(Arc::clone(&slot.modification)),
                    position => {
                        residues.insert(usize::from(position), Arc::clone(&slot.modification));
                    }
                }
            }
        }
        let residues: Vec<(usize, Arc<ResidueModification>)> = residues.into_iter().collect();
        seq.with_resolved_modifications(&residues, n_term, c_term)
    }

    /// The realized length of an SNES mother whose (M+H)+ is closest to
    /// `target_mh_plus` within `[-lower, +upper]` (source
    /// `realizeSNESLength`), or `None` (the source's `-1`) when no length in
    /// `[peptide:min_size, mother length]` fits, or outside SNES mode.
    ///
    /// Residues are accumulated from the anchored terminus (left to right for
    /// Single-N, right to left for Single-C) with fixed modification deltas,
    /// on top of water, a proton and both fixed terminal deltas. Ppm
    /// magnitudes are converted at `max(target, 0)`. Ties keep the shorter
    /// length.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for coordinates outside
    /// `fasta_entries`.
    pub fn realize_snes_length(
        &self,
        mother: &Peptide,
        fasta_entries: &[FASTAEntry],
        target_mh_plus: f64,
        tolerance_lower_magnitude: f64,
        tolerance_upper_magnitude: f64,
        tolerance_ppm: bool,
    ) -> Result<Option<usize>> {
        if !self.is_snes_mode() {
            return Ok(None);
        }
        let (protein, start, len) =
            coordinates(mother, fasta_entries, usize::from(mother.sequence.1))?;
        let masses = mass_tables()?;
        let tables = &self.tables;
        let single_c = Self::is_single_c_mother(mother.mod_bitmask);
        let base =
            masses.water + PROTON_MASS_U + tables.fixed_nterm_delta + tables.fixed_cterm_delta;
        let (lo, hi) = if tolerance_ppm {
            (
                tolerance_lower_magnitude * target_mh_plus.max(0.0) * 1e-6,
                tolerance_upper_magnitude * target_mh_plus.max(0.0) * 1e-6,
            )
        } else {
            (tolerance_lower_magnitude, tolerance_upper_magnitude)
        };
        let bytes = protein.as_bytes();
        let mut cumulative = 0.0;
        let mut best_abs_delta = f64::MAX;
        let mut best = None;
        for k in 0..len {
            let index = if single_c {
                start + (len - 1 - k)
            } else {
                start + k
            };
            let aa = usize::from(bytes[index]);
            cumulative += masses.residue[aa] + tables.fixed_deltas[aa];
            let length = k + 1;
            if length < self.params.peptide_min_length {
                continue;
            }
            let delta = base + cumulative - target_mh_plus;
            if delta < -lo || delta > hi {
                continue;
            }
            if delta.abs() < best_abs_delta {
                best_abs_delta = delta.abs();
                best = Some(length);
            }
        }
        Ok(best)
    }

    /// A realized SNES sub-peptide as a sequence (source
    /// `reconstructRealizedSubSequence`): the first `realized_length`
    /// residues of a Single-N mother or the last ones of a Single-C mother,
    /// with fixed modifications and the variable modifications of
    /// `subset_bitmask` (slots of the realized sub-peptide; 0 = none).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `realized_length` exceeds the
    /// mother or coordinates fall outside `fasta_entries`.
    pub fn reconstruct_realized_sub_sequence(
        &self,
        mother: &Peptide,
        fasta_entries: &[FASTAEntry],
        realized_length: usize,
        subset_bitmask: u32,
    ) -> Result<AASequence> {
        let mother_length = usize::from(mother.sequence.1);
        if realized_length > mother_length {
            return Err(Error::InvalidValue(
                "realized length exceeds the mother".into(),
            ));
        }
        let (protein, mother_start, _) = coordinates(mother, fasta_entries, mother_length)?;
        let start = if Self::is_single_c_mother(mother.mod_bitmask) {
            mother_start + mother_length - realized_length
        } else {
            mother_start
        };
        self.sequence_with_mods(protein, start, realized_length, subset_bitmask, false)
    }

    /// The half-open peptide index range whose precursor m/z lies in
    /// `[precursor_mass + window.0, precursor_mass + window.1]` (source
    /// `getPeptidesInMassWindow`), computed in `f32`. A reversed window gives
    /// `(0, 0)`.
    pub fn get_peptides_in_mass_window(
        &self,
        precursor_mass: f32,
        window: (f32, f32),
    ) -> (usize, usize) {
        if window.0 > window.1 {
            return (0, 0);
        }
        let low = precursor_mass + window.0;
        let high = precursor_mass + window.1;
        let left = self.peptides.partition_point(|p| p.precursor_mz < low);
        let right = self.peptides.partition_point(|p| !(high < p.precursor_mz));
        (left, right)
    }

    /// The signed window `(lo <= 0, hi >= 0)` around `precursor_mass`, in Da
    /// (source private `computeMassWindow_`); ppm magnitudes are converted in
    /// `f32` at `precursor_mass`.
    pub fn compute_mass_window(&self, precursor_mass: f32) -> (f32, f32) {
        let p = &self.params;
        if p.precursor_mass_tolerance_unit_ppm {
            (
                -ppm_to_mass_f32(p.precursor_mass_tolerance_lower as f32, precursor_mass),
                ppm_to_mass_f32(p.precursor_mass_tolerance_upper as f32, precursor_mass),
            )
        } else {
            (
                -(p.precursor_mass_tolerance_lower as f32),
                p.precursor_mass_tolerance_upper as f32,
            )
        }
    }

    /// Bucket range to walk for a mass and tolerance: `lower_bound` of
    /// `mass - tol`, one bucket back, to `upper_bound` of `mass + tol`.
    fn bucket_range(&self, low: f32, high: f32) -> (usize, usize) {
        let mut left = self.bucket_min_mz.partition_point(|&m| m < low);
        let right = self.bucket_min_mz.partition_point(|&m| m <= high);
        left = left.saturating_sub(1);
        (left, right.max(left))
    }

    fn bucket(&self, j: usize) -> &[Fragment] {
        let start = j * self.bucket_size;
        let end = ((j + 1) * self.bucket_size).min(self.fragments.len());
        &self.fragments[start.min(end)..end]
    }

    /// Peak mass at a charge and the fragment tolerance in Da, with the
    /// source's float semantics.
    fn adjusted_mass_and_tolerance(&self, peak: &Peak1D, charge: u16) -> (f32, f32) {
        let adjusted = (peak.mz * f64::from(charge as f32)
            - f64::from(i32::from(charge) - 1) * PROTON_MASS_U) as f32;
        let p = &self.params;
        let tolerance = if p.fragment_mz_tolerance_unit_ppm {
            ppm_to_mass_f32(p.fragment_mz_tolerance, adjusted)
        } else {
            p.fragment_mz_tolerance
        };
        (adjusted, tolerance)
    }

    /// Fragments matching one peak within the peptide range `[first,
    /// second)` (source `query`): the peak is converted to its charge-1 mass,
    /// `mz * charge - (charge - 1) * PROTON_MASS_U` in `f32`, and matched
    /// within the fragment tolerance (converted at that mass for ppm).
    pub fn query(
        &self,
        peak: &Peak1D,
        peptide_idx_range: (usize, usize),
        peak_charge: u16,
    ) -> Vec<Hit> {
        let (adjusted, tolerance) = self.adjusted_mass_and_tolerance(peak, peak_charge);
        let (begin, end) = self.bucket_range(adjusted - tolerance, adjusted + tolerance);
        let mut hits = Vec::new();
        for j in begin..end {
            let slice = self.bucket(j);
            let from = slice.partition_point(|f| (f.peptide_idx as usize) < peptide_idx_range.0);
            for f in &slice[from..] {
                if f.peptide_idx as usize >= peptide_idx_range.1 {
                    break;
                }
                if adjusted >= f.fragment_mz - tolerance && adjusted <= f.fragment_mz + tolerance {
                    hits.push(Hit {
                        peptide_idx: f.peptide_idx,
                        fragment_mz: f.fragment_mz,
                    });
                }
            }
        }
        hits
    }

    /// Count fragment matches for one (precursor charge, isotope error)
    /// block and append the candidates reaching `fragment:min_matched_ions`
    /// (at least 1) in ascending peptide index (source private
    /// `queryPeaks`). Peaks are walked at fragment charges `1..=min(charge,
    /// fragment:max_charge)`.
    fn query_peaks(
        &self,
        candidates: &mut SpectrumMatchesTopN,
        spectrum: &MSSpectrum,
        range: (usize, usize),
        isotope_error: i16,
        precursor_charge: u16,
    ) {
        if range.0 >= range.1 {
            return;
        }
        let actual_max = precursor_charge.min(self.params.max_fragment_charge);
        let mut counts: BTreeMap<usize, u32> = BTreeMap::new();
        for peak in &spectrum.peaks {
            for fragment_charge in 1..=actual_max {
                let (adjusted, tolerance) = self.adjusted_mass_and_tolerance(peak, fragment_charge);
                let (begin, end) = self.bucket_range(adjusted - tolerance, adjusted + tolerance);
                for j in begin..end {
                    let slice = self.bucket(j);
                    let from = slice.partition_point(|f| (f.peptide_idx as usize) < range.0);
                    for f in &slice[from..] {
                        if f.peptide_idx as usize >= range.1 {
                            break;
                        }
                        if adjusted >= f.fragment_mz - tolerance
                            && adjusted <= f.fragment_mz + tolerance
                        {
                            let cell = counts.entry(f.peptide_idx as usize).or_insert(0);
                            *cell = cell.saturating_add(1);
                        }
                    }
                }
            }
        }
        let emit_min = u32::from(self.params.min_matched_peaks).max(1);
        for (peptide_idx, count) in counts {
            if count >= emit_min {
                candidates.hits.push(SpectrumMatch {
                    num_matched: count,
                    precursor_charge,
                    isotope_error,
                    peptide_idx,
                    ..SpectrumMatch::default()
                });
            }
        }
    }

    /// Rank candidates and keep the best `scoring:max_candidates_per_spectrum`
    /// with at least `fragment:min_matched_ions` matches (source private
    /// `trimHits`). Ranking: more matches, smaller |isotope error|, smaller
    /// isotope error, smaller charge, smaller peptide index.
    pub fn trim_hits(&self, hits: &mut SpectrumMatchesTopN) {
        hits.hits.sort_by(by_rank);
        hits.hits.truncate(self.params.max_processed_hits as usize);
        let min = u32::from(self.params.min_matched_peaks);
        if hits.hits.first().is_some_and(|h| h.num_matched < min) {
            hits.hits.clear();
        }
        if let Some(last_kept) = hits.hits.iter().rposition(|h| h.num_matched >= min) {
            hits.hits.truncate(last_kept + 1);
        }
    }

    /// Query every isotope error of a closed search (or only 0 in open
    /// search) for one precursor charge (source private
    /// `searchDifferentPrecursorRanges`).
    fn search_different_precursor_ranges(
        &self,
        spectrum: &MSSpectrum,
        precursor_mass: f32,
        sms: &mut SpectrumMatchesTopN,
        charge: u16,
    ) {
        let p = &self.params;
        let open = p.is_open_search_mode();
        let (lo, hi) = if open {
            (0, 0)
        } else {
            (p.min_isotope_error, p.max_isotope_error)
        };
        for isotope_error in lo..=hi {
            let shifted = precursor_mass + f32::from(isotope_error) * (C13C12_MASSDIFF_U as f32);
            let window = self.compute_mass_window(shifted);
            let range = self.get_peptides_in_mass_window(shifted, window);
            self.query_peaks(sms, spectrum, range, isotope_error, charge);
        }
    }

    /// Query one MS2 spectrum against the index (source `querySpectrum`),
    /// appending the trimmed candidates to `sms`.
    ///
    /// The precursor charge is taken from the single precursor, or every
    /// charge of `precursor:min_charge..=max_charge` when it is 0; each
    /// charge's (M+H)+ is `mz * charge - (charge - 1) * PROTON_MASS_U` in
    /// `f32`. Candidates of all charges and isotope errors are then ranked
    /// together by [`FragmentIndex::trim_hits`]. Nothing is appended for an
    /// unbuilt index (with a logged warning), an empty or non-MS2 spectrum,
    /// or a spectrum without exactly one precursor (warning).
    ///
    /// # Errors
    ///
    /// Returns [`Error::Unsupported`] in SNES mode: the source's
    /// `querySpectrumSNES_` (byte-count pre-pass and fragment-bin precursor
    /// filter) is not ported yet.
    pub fn query_spectrum(
        &self,
        spectrum: &MSSpectrum,
        sms: &mut SpectrumMatchesTopN,
    ) -> Result<()> {
        if !self.is_build {
            let _ = crate::openms_log_warn!("FragmentIndex not yet build ");
            return Ok(());
        }
        if spectrum.peaks.is_empty() || spectrum.ms_level != 2 {
            return Ok(());
        }
        if spectrum.precursors.len() != 1 {
            let _ = crate::openms_log_warn!("Number of precursors is not equal 1 ");
            return Ok(());
        }
        if self.is_snes_mode() {
            return Err(Error::Unsupported(
                "FragmentIndex SNES spectrum query (querySpectrumSNES_) is not ported".into(),
            ));
        }
        let precursor = &spectrum.precursors[0];
        let charges: Vec<u16> = if precursor.charge != 0 {
            // The source pushes the int charge into a size_t vector and
            // narrows it to uint16_t, so a negative charge wraps.
            vec![precursor.charge as u16]
        } else {
            (self.params.min_precursor_charge..=self.params.max_precursor_charge).collect()
        };
        for charge in charges {
            let scaled = (precursor.mz as f32) * f32::from(charge);
            let mz = (f64::from(scaled) - f64::from(i32::from(charge) - 1) * PROTON_MASS_U) as f32;
            self.search_different_precursor_ranges(spectrum, mz, sms, charge);
        }
        self.trim_hits(sms);
        Ok(())
    }

    /// As [`FragmentIndex::query_spectrum`], with the FASTA entries the SNES
    /// query with variable modifications needs (source three-argument
    /// `querySpectrum`). Outside SNES mode the entries are ignored.
    ///
    /// # Errors
    ///
    /// As [`FragmentIndex::query_spectrum`].
    pub fn query_spectrum_with_fasta(
        &self,
        spectrum: &MSSpectrum,
        fasta_entries: &[FASTAEntry],
        sms: &mut SpectrumMatchesTopN,
    ) -> Result<()> {
        let _ = fasta_entries;
        self.query_spectrum(spectrum, sms)
    }
}

fn slot(position: u16, e: &VarModEntry) -> ModSlot {
    ModSlot {
        position,
        delta_mass: e.delta_mass,
        modification: Arc::clone(&e.modification),
    }
}

/// Source `Math::ppmToMass<float>`: `(ppm / 1e6) * mz_ref` in `f32`.
fn ppm_to_mass_f32(ppm: f32, mz_ref: f32) -> f32 {
    (ppm / 1e6_f32) * mz_ref
}

fn by_rank(a: &SpectrumMatch, b: &SpectrumMatch) -> Ordering {
    b.num_matched
        .cmp(&a.num_matched)
        .then_with(|| {
            a.isotope_error
                .unsigned_abs()
                .cmp(&b.isotope_error.unsigned_abs())
        })
        .then_with(|| a.isotope_error.cmp(&b.isotope_error))
        .then_with(|| a.precursor_charge.cmp(&b.precursor_charge))
        .then_with(|| a.peptide_idx.cmp(&b.peptide_idx))
}

/// Stable sort by `(precursor_mz, protein_idx)`.
fn sort_peptides(peptides: &mut [Peptide]) {
    peptides.sort_by(|a, b| {
        a.precursor_mz
            .total_cmp(&b.precursor_mz)
            .then(a.protein_idx.cmp(&b.protein_idx))
    });
}

fn coordinates<'a>(
    peptide: &Peptide,
    fasta_entries: &'a [FASTAEntry],
    length: usize,
) -> Result<(&'a str, usize, usize)> {
    let protein = fasta_entries
        .get(peptide.protein_idx as usize)
        .ok_or_else(|| Error::InvalidValue("peptide protein index out of range".into()))?;
    let start = usize::from(peptide.sequence.0);
    if start + length > protein.sequence.len() || !protein.sequence.is_char_boundary(start) {
        return Err(Error::InvalidValue(
            "peptide coordinates out of range".into(),
        ));
    }
    Ok((&protein.sequence, start, length))
}

/// Source `initModificationTables_`.
fn init_modification_tables(p: &FragmentIndexParameters) -> Result<ModTables> {
    let mut tables = ModTables {
        fixed_deltas: [0.0; 128],
        fixed_mods: vec![None; 128],
        fixed_nterm_delta: 0.0,
        fixed_cterm_delta: 0.0,
        fixed_nterm_mod: None,
        fixed_cterm_mod: None,
        variable: vec![Vec::new(); 128],
        variable_nterm: Vec::new(),
        variable_cterm: Vec::new(),
    };
    let resolve = |names: &[String]| -> Result<Vec<Arc<ResidueModification>>> {
        let names: Vec<&str> = names.iter().map(String::as_str).collect();
        ModifiedPeptideGenerator::get_modifications(&names)
    };
    let terminal_only =
        |m: &ResidueModification| matches!(m.origin(), None | Some('X') | Some('.'));
    let is_n =
        |t: TermSpecificity| matches!(t, TermSpecificity::NTerm | TermSpecificity::ProteinNTerm);
    let is_c =
        |t: TermSpecificity| matches!(t, TermSpecificity::CTerm | TermSpecificity::ProteinCTerm);
    for m in resolve(&p.modifications_fixed)? {
        let term = m.term_specificity();
        let delta = m.diff_mono_mass();
        if !terminal_only(&m) && term == TermSpecificity::Anywhere {
            let origin = m.origin().map_or(0, |c| c as usize);
            if origin < 128 {
                tables.fixed_deltas[origin] = delta;
                tables.fixed_mods[origin] = Some(Arc::clone(&m));
            }
        } else if is_n(term) {
            tables.fixed_nterm_delta = delta;
            tables.fixed_nterm_mod = Some(Arc::clone(&m));
        } else if is_c(term) {
            tables.fixed_cterm_delta = delta;
            tables.fixed_cterm_mod = Some(Arc::clone(&m));
        }
    }
    for m in resolve(&p.modifications_variable)? {
        let term = m.term_specificity();
        let entry = VarModEntry {
            delta_mass: m.diff_mono_mass(),
            modification: Arc::clone(&m),
            term_spec: term,
        };
        if terminal_only(&m) {
            if is_n(term) {
                tables.variable_nterm.push(entry);
            } else if is_c(term) {
                tables.variable_cterm.push(entry);
            }
        } else {
            let origin = m.origin().map_or(0, |c| c as usize);
            if origin < 128 {
                tables.variable[origin].push(entry);
            }
        }
    }
    Ok(tables)
}
