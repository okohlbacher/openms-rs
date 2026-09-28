// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Subset-neighbor peptide search, ported from `ANALYSIS/ID/NeighborSeq.h`.
//!
//! Find peptides from a wider pool (typically a FASTA digest) that are
//! spectral neighbors of a smaller "relevant" peptide set, useful when only
//! part of a complex sample is of interest. Two peptides are neighbors when
//! their precursor masses are within tolerance and their theoretical b/y
//! fragment spectra share enough peaks.
//!
//! Background: Cormen et al., *J. Proteome Research* 2021,
//! `10.1021/acs.jproteome.1c00483`
//! (<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8489664/>).
//!
//! See `docs/NEIGHBOR_SEQ_SUPPORT.md` for the API mapping and evidence.

use crate::chemistry::AASequence;
use crate::chemistry::theoretical::{TheoreticalIonSeries, TheoreticalSpectrumGenerator};
use crate::concept::math_functions::ppm_to_mass;
use crate::kernel::MSSpectrum;
use crate::{Error, Result};

/// Subset-neighbor peptide search (source `NeighborSeq`).
///
/// Constructed once with the relevant peptides, then queried with
/// [`NeighborSeq::is_neighbor_peptide`] for each candidate. After the queries,
/// [`NeighborSeq::neighbor_stats`] summarises how many relevant peptides had
/// zero, one, or multiple neighbors.
///
/// ```
/// use openms::analysis::id::neighbor_seq::NeighborSeq;
/// use openms::chemistry::AASequence;
/// let mut ns = NeighborSeq::new(vec![AASequence::parse("VGEFK")?])?;
/// assert!(ns.is_neighbor_peptide(&AASequence::parse("GLDFK")?, 0.01, false, 0.25, 0.05)?);
/// assert_eq!(ns.neighbor_stats().findable_one_neighbor, 1);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug)]
pub struct NeighborSeq {
    /// Digested relevant peptides, owned (source `digested_relevant_peptides_`).
    digested_relevant_peptides: Vec<AASequence>,
    /// Distinct masses in ascending order, each with the positions in
    /// `digested_relevant_peptides` that have exactly this mass, in insertion
    /// order (source `mass_position_map_`, a `std::map<double, vector<int>>`).
    mass_positions: Vec<(f64, Vec<usize>)>,
    /// Charge-1 b/y generator including the b1 prefix ion (source `spec_gen_`).
    spec_gen: TheoreticalSpectrumGenerator,
    /// Neighbors found per relevant peptide; `None` marks a peptide containing
    /// `X`, which cannot be searched (the source stores `-1`).
    neighbor_stats: Vec<Option<usize>>,
}

impl NeighborSeq {
    /// Construct from the "relevant" digested peptides.
    ///
    /// Builds the internal mass index and configures the theoretical-spectrum
    /// generator for charge-1 b/y ions including the b1 prefix ion. Relevant
    /// peptides containing the unknown-amino-acid residue `X` are skipped and
    /// counted as unfindable; the count is logged as a warning, as the source
    /// does (always, even when it is zero).
    ///
    /// The source takes the vector by rvalue reference and owns it for its
    /// lifetime; this takes it by value, which is the same ownership contract.
    ///
    /// # Errors
    ///
    /// Returns an error when a relevant peptide without `X` has no known
    /// monoisotopic mass (for example the ambiguous residues `B`/`Z`, or an
    /// annotation without a mass). The source computes a mass from its residue
    /// table for those; this port refuses rather than guess.
    pub fn new(digested_relevant_peptides: Vec<AASequence>) -> Result<Self> {
        let spec_gen = TheoreticalSpectrumGenerator {
            ion_series: vec![TheoreticalIonSeries::B, TheoreticalIonSeries::Y],
            add_first_prefix_ion: true,
            ..TheoreticalSpectrumGenerator::default()
        };
        let mut neighbor_stats = vec![Some(0); digested_relevant_peptides.len()];
        let mut by_mass: Vec<(f64, usize)> = Vec::with_capacity(digested_relevant_peptides.len());
        let mut skipped = 0_usize;
        for (index, peptide) in digested_relevant_peptides.iter().enumerate() {
            if contains_unknown_residue(peptide) {
                neighbor_stats[index] = None;
                skipped += 1;
                continue;
            }
            by_mass.push((peptide.mono_mass()?, index));
        }
        // A stable sort keeps insertion order within one mass, as
        // `map[mass].push_back(i)` does.
        by_mass.sort_by(|left, right| left.0.total_cmp(&right.0));
        let mut mass_positions: Vec<(f64, Vec<usize>)> = Vec::new();
        for (mass, index) in by_mass {
            match mass_positions.last_mut() {
                Some((last, positions)) if *last == mass => positions.push(index),
                _ => mass_positions.push((mass, vec![index])),
            }
        }
        let _ = crate::openms_log_warn!(
            "Skipped {skipped}/{} peptides with unknown('X') amino acids.",
            digested_relevant_peptides.len()
        );
        Ok(Self {
            digested_relevant_peptides,
            mass_positions,
            spec_gen,
            neighbor_stats,
        })
    }

    /// The relevant peptides this instance was constructed with.
    pub fn relevant_peptides(&self) -> &[AASequence] {
        &self.digested_relevant_peptides
    }

    /// Theoretical spectrum of `peptide_sequence` with b/y ions at charge 1.
    ///
    /// Includes all b and y ions with charge 1, even the prefix ion b1, but no
    /// internal ions; peaks are sorted by m/z.
    ///
    /// # Errors
    ///
    /// As [`TheoreticalSpectrumGenerator::generate`], notably when a residue
    /// has no known mass.
    pub fn generate_spectrum(&self, peptide_sequence: &AASequence) -> Result<MSSpectrum> {
        self.spec_gen.generate(peptide_sequence, 1, 1, None)
    }

    /// Whether two spectra share enough peaks, in `mz_bin_size` m/z bins, to
    /// be considered neighbors.
    ///
    /// Computes the Dice-style fraction `2 * shared_peaks / (spec1.len() +
    /// spec2.len())` and returns `true` when it is **strictly** greater than
    /// `min_shared_ion_fraction` (in `[0, 1]`). All peaks of both spectra are
    /// considered; shared peaks come from [`NeighborSeq::compute_shared_ion_count`].
    /// Typical bin sizes are `0.05` Th for high-resolution and `1.0005079` Th
    /// for low-resolution data. Two empty spectra give `0/0`, which is not
    /// greater than any threshold, so the result is `false`, as in the source.
    pub fn is_neighbor_spectrum(
        spec1: &MSSpectrum,
        spec2: &MSSpectrum,
        min_shared_ion_fraction: f64,
        mz_bin_size: f64,
    ) -> bool {
        let shared = Self::compute_shared_ion_count(spec1, spec2, mz_bin_size);
        let fraction_shared =
            (2.0 * shared as f64) / (spec1.peaks.len() + spec2.peaks.len()) as f64;
        fraction_shared > min_shared_ion_fraction
    }

    /// Number of shared ions between two spectra.
    ///
    /// Every peak is mapped to the bin `trunc(mz / mz_bin_size)` (bins are
    /// `[a, b[` for positive m/z), and the two bin sequences are intersected
    /// with a sorted-range merge: a matched peak of `spec1` is consumed, and
    /// each peak of `spec2` is consumed once visited. All peaks are
    /// considered. Use [`NeighborSeq::generate_spectrum`] to obtain theoretical
    /// b/y spectra.
    ///
    /// Both spectra must be sorted by m/z; unsorted input gives a meaningless
    /// count (the merge always terminates). The source truncates to a 32-bit
    /// `int`, which is undefined for bins beyond its range; this truncates to
    /// `i64` with Rust's saturating float conversion, identical within the
    /// source's defined range.
    pub fn compute_shared_ion_count(
        spec1: &MSSpectrum,
        spec2: &MSSpectrum,
        mz_bin_size: f64,
    ) -> usize {
        let bin = |mz: f64| (mz / mz_bin_size) as i64;
        let (mut first1, mut first2) =
            (spec1.peaks.iter().peekable(), spec2.peaks.iter().peekable());
        let mut count = 0;
        while let (Some(p1), Some(p2)) = (first1.peek(), first2.peek()) {
            let (val1, val2) = (bin(p1.mz), bin(p2.mz));
            if val1 < val2 {
                first1.next();
            } else {
                if val1 == val2 {
                    first1.next();
                    count += 1;
                }
                first2.next();
            }
        }
        count
    }

    /// Whether `neighbor_candidate` is a spectral neighbor of any relevant
    /// peptide.
    ///
    /// Looks up the relevant peptides whose precursor mass is within
    /// `mass_tolerance_pc` of the candidate's monoisotopic mass (inclusive at
    /// both ends), generates b/y spectra for each of them and the candidate,
    /// and compares them with [`NeighborSeq::is_neighbor_spectrum`] using
    /// `min_shared_ion_fraction` and `mz_bin_size`. Returns `true` when **any**
    /// relevant peptide qualifies, but keeps iterating so the per-peptide
    /// neighbor counters are updated for every match. Call
    /// [`NeighborSeq::neighbor_stats`] once all candidates have been queried.
    ///
    /// # Arguments
    ///
    /// * `mass_tolerance_pc` - maximum precursor mass difference, in Da, or in
    ///   ppm when `mass_tolerance_pc_ppm` is `true`.
    /// * `mz_bin_size` - typically `0.05` Th (high-res) or `1.0005079` Th
    ///   (low-res).
    ///
    /// The ppm conversion calls the source's `Math::ppmToMass(mono_weight,
    /// tolerance)` with its two arguments swapped relative to the declared
    /// `(ppm, mz_ref)`; the product is the same up to rounding, and this port
    /// keeps the source's operand order so the window is bit-identical.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `mass_tolerance_pc` is negative or
    /// not finite (the source only `assert`s non-negativity in debug builds),
    /// and propagates mass or spectrum-generation errors, for example for a
    /// candidate containing a residue without a known mass. Counters are
    /// updated only after every comparison succeeded.
    pub fn is_neighbor_peptide(
        &mut self,
        neighbor_candidate: &AASequence,
        mass_tolerance_pc: f64,
        mass_tolerance_pc_ppm: bool,
        min_shared_ion_fraction: f64,
        mz_bin_size: f64,
    ) -> Result<bool> {
        let (from, to) = self.find_candidate_positions(
            neighbor_candidate.mono_mass()?,
            mass_tolerance_pc,
            mass_tolerance_pc_ppm,
        )?;
        if from == to {
            return Ok(false);
        }
        let spec = self.generate_spectrum(neighbor_candidate)?;
        let mut hits = Vec::new();
        for (_, positions) in &self.mass_positions[from..to] {
            for &index in positions {
                let neighbor_spec =
                    self.generate_spectrum(&self.digested_relevant_peptides[index])?;
                if Self::is_neighbor_spectrum(
                    &spec,
                    &neighbor_spec,
                    min_shared_ion_fraction,
                    mz_bin_size,
                ) {
                    hits.push(index);
                }
            }
        }
        for &index in &hits {
            if let Some(count) = &mut self.neighbor_stats[index] {
                *count += 1;
            }
        }
        Ok(!hits.is_empty())
    }

    /// Statistics of how many neighbors were found per relevant peptide by the
    /// [`NeighborSeq::is_neighbor_peptide`] calls so far (source
    /// `getNeighborStats`).
    pub fn neighbor_stats(&self) -> NeighborStats {
        let mut stats = NeighborStats::default();
        for count in &self.neighbor_stats {
            match count {
                None => stats.unfindable_peptides += 1,
                Some(0) => stats.findable_no_neighbors += 1,
                Some(1) => stats.findable_one_neighbor += 1,
                Some(_) => stats.findable_multiple_neighbors += 1,
            }
        }
        stats
    }

    /// Index range `[from, to)` into `mass_positions` whose masses lie within
    /// the tolerance window around `mono_weight` (source
    /// `findCandidatePositions_`, a `lower_bound`/`upper_bound` pair).
    fn find_candidate_positions(
        &self,
        mono_weight: f64,
        mass_tolerance: f64,
        mass_tolerance_pc_ppm: bool,
    ) -> Result<(usize, usize)> {
        if !(mass_tolerance.is_finite() && mass_tolerance >= 0.0) {
            return Err(Error::InvalidValue(
                "precursor mass tolerance must be finite and non-negative".into(),
            ));
        }
        let tolerance = if mass_tolerance_pc_ppm {
            // Source argument order, see `is_neighbor_peptide`.
            ppm_to_mass(mono_weight, mass_tolerance)?
        } else {
            mass_tolerance
        };
        let (low, high) = (mono_weight - tolerance, mono_weight + tolerance);
        let from = self.mass_positions.partition_point(|(mass, _)| *mass < low);
        let to = self
            .mass_positions
            .partition_point(|(mass, _)| *mass <= high);
        Ok((from, to.max(from)))
    }
}

/// Whether the sequence contains the unknown residue `X` (source
/// `AASequence::has(*ResidueDB::getResidue('X'))`).
fn contains_unknown_residue(peptide: &AASequence) -> bool {
    peptide.as_str().contains('X')
}

/// Statistics of how many neighbors were found per relevant peptide (source
/// `NeighborSeq::NeighborStats`).
///
/// The four counters are mutually exclusive categories.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NeighborStats {
    /// Relevant peptides containing an `X` (unknown amino acid), which cannot
    /// be searched for neighbors.
    pub unfindable_peptides: usize,
    /// Relevant peptides with no neighbors.
    pub findable_no_neighbors: usize,
    /// Relevant peptides with exactly one neighbor.
    pub findable_one_neighbor: usize,
    /// Relevant peptides with multiple neighbors.
    pub findable_multiple_neighbors: usize,
}

impl NeighborStats {
    /// Sum of all four categories, i.e. the number of relevant peptides
    /// registered at construction time. Saturates instead of overflowing.
    pub fn total(&self) -> usize {
        self.unfindable_peptides
            .saturating_add(self.findable_no_neighbors)
            .saturating_add(self.findable_one_neighbor)
            .saturating_add(self.findable_multiple_neighbors)
    }

    /// Percentage (0..100, truncated) of `count` relative to
    /// [`NeighborStats::total`]; 0 when the total is 0, avoiding a division by
    /// zero (source `percentOfTotal_`). Computed in 128-bit arithmetic, so
    /// `count * 100` cannot overflow as the source's `int` product can.
    pub fn percent_of_total(&self, count: usize) -> u128 {
        match self.total() {
            0 => 0,
            total => (count as u128 * 100) / total as u128,
        }
    }

    fn format(&self, count: usize) -> String {
        format!("{count} ({}%)", self.percent_of_total(count))
    }

    /// [`NeighborStats::unfindable_peptides`] formatted as `"X (Y%)"`; `"X (0%)"`
    /// when the total is 0.
    pub fn unfindable(&self) -> String {
        self.format(self.unfindable_peptides)
    }

    /// [`NeighborStats::findable_no_neighbors`] formatted as `"X (Y%)"` (source
    /// `noNB`); `"X (0%)"` when the total is 0.
    pub fn no_nb(&self) -> String {
        self.format(self.findable_no_neighbors)
    }

    /// [`NeighborStats::findable_one_neighbor`] formatted as `"X (Y%)"` (source
    /// `oneNB`); `"X (0%)"` when the total is 0.
    pub fn one_nb(&self) -> String {
        self.format(self.findable_one_neighbor)
    }

    /// [`NeighborStats::findable_multiple_neighbors`] formatted as `"X (Y%)"`
    /// (source `multiNB`); `"X (0%)"` when the total is 0.
    pub fn multi_nb(&self) -> String {
        self.format(self.findable_multiple_neighbors)
    }
}
