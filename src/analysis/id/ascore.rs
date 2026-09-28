// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/AScore.h`: phosphorylation-site localisation with the
//! AScore of Beausoleil et al. (Nat. Biotechnol. 2006).
//!
//! For a peptide hit and its MS/MS spectrum, [`AScore::compute`] enumerates
//! every placement of the hit's phosphorylations on its S/T/Y residues, scores
//! each placement against the spectrum's top peaks per 100 Th window, and
//! reports for each site of the best placement the score difference to the
//! best placement without that site, from the ions that tell the two apart.
//!
//! The source's protected helpers are public here, because they are the
//! surface its class test drives and the steps a reader needs to follow the
//! algorithm; each is documented as a step of [`AScore::compute`]. The source's
//! mutable `base_match_probability_` member is replaced by an explicit
//! argument, so every method takes `&self`. The API mapping, preserved
//! conventions and native differences are in `docs/PSCORE_ASCORE_SUPPORT.md`.

use crate::chemistry::AASequence;
use crate::chemistry::theoretical::TheoreticalSpectrumGenerator;
use crate::comparison::{Tolerance, matched_alignment};
use crate::concept::math_functions::binomial_cdf_complement;
use crate::identification::PeptideHit;
use crate::metadata::MetaValue;
use crate::param::{DefaultParamHandler, Param, ParamEntry, ParamNode, ParamValue};
use crate::{Error, MSSpectrum, Peak1D, Result};
use std::cmp::Ordering;

/// Source `DefaultParamHandler` name of the AScore.
pub const ASCORE_NAME: &str = "AScore";

/// Largest number of phosphosite permutations [`AScore::compute`] enumerates.
///
/// The source enumerates every placement before comparing the count with
/// `max_num_perm`, and `max_num_perm = 0` means "no restriction", so a long
/// S/T/Y-rich peptide can ask for billions of placements. This port counts the
/// placements first, returns the source's aborted result when the count exceeds
/// a positive `max_num_perm`, and refuses a count above this ceiling whatever
/// `max_num_perm` says. The default `max_num_perm` of 16384 lies well below it.
pub const MAX_ASCORE_PERMUTATIONS: usize = 1 << 20;

/// Largest number of 100 Th windows [`AScore::peak_picking_per_windows_in_spectrum`]
/// allocates (a spectrum spanning 10^7 Th). The source has no ceiling.
pub const MAX_ASCORE_WINDOWS: usize = 100_000;

/// The number of peak depths (1 to 10 peaks per window) every placement is
/// scored at.
pub const ASCORE_PEAK_DEPTHS: usize = 10;

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

fn count(value: usize) -> Result<u32> {
    u32::try_from(value).map_err(|_| bad("AScore ion count exceeds u32::MAX"))
}

/// One phosphorylation site of the best placement and its runner-up
/// (source `ProbablePhosphoSites`).
///
/// Every field is an index; the source value-initialises them to zero, which
/// `Default` reproduces.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ProbablePhosphoSites {
    /// Residue index of this site in the best placement.
    pub first: usize,
    /// Residue index of the site the runner-up placement uses instead.
    pub second: usize,
    /// Index of the best placement, with this site phosphorylated.
    pub seq_1: usize,
    /// Index of the best placement with this site unphosphorylated.
    pub seq_2: usize,
    /// The peak depth (filtering level, 1 to 10) that gave the maximum
    /// discriminatory score.
    pub peak_depth: usize,
    /// Source member `AScore`. The source never assigns it; it stays zero and
    /// the per-site scores are reported as `AScore_<n>` metadata instead.
    pub a_score: usize,
}

/// The AScore phosphosite localisation (source `AScore`, a
/// `DefaultParamHandler` named [`ASCORE_NAME`]).
///
/// The fields are the source members set by `updateMembers_`; build one with
/// [`AScore::from_param`] to get the source's parameter checks, or set the
/// fields directly. [`AScore::compute`] re-checks what it relies on.
#[derive(Clone, Debug, PartialEq)]
pub struct AScore {
    /// Fragment mass tolerance for spectrum comparisons (`fragment_mass_tolerance`,
    /// default 0.05, at least 0).
    pub fragment_mass_tolerance: f64,
    /// Whether the tolerance is in ppm rather than Da (`fragment_mass_unit`,
    /// `Da` or `ppm`, default `Da`).
    pub fragment_tolerance_ppm: bool,
    /// Longest peptide, in residues, that is scored; 0 means no restriction
    /// (`max_peptide_length`, default 40, advanced).
    pub max_peptide_length: usize,
    /// Largest number of placements that is scored; 0 means no restriction
    /// (`max_num_perm`, default 16384, advanced). [`MAX_ASCORE_PERMUTATIONS`]
    /// applies regardless.
    pub max_permutations: usize,
    /// Score of an unambiguous assignment, where every S/T/Y is phosphorylated
    /// (`unambiguous_score`, default 1000, advanced). An integer parameter in
    /// the source, converted to its `double` member.
    pub unambiguous_score: i32,
}

impl Default for AScore {
    fn default() -> Self {
        Self {
            fragment_mass_tolerance: 0.05,
            fragment_tolerance_ppm: false,
            max_peptide_length: 40,
            max_permutations: 16384,
            unambiguous_score: 1000,
        }
    }
}

fn entry(name: &str, value: ParamValue, description: &str, advanced: bool) -> ParamEntry {
    let mut entry = ParamEntry {
        name: name.into(),
        description: description.into(),
        value,
        ..ParamEntry::default()
    };
    if advanced {
        entry.tags.insert("advanced".into());
    }
    entry
}

fn source_int(value: usize, key: &str) -> Result<i64> {
    i32::try_from(value)
        .map(i64::from)
        .map_err(|_| bad(&format!("'{key}' does not fit the source 32-bit parameter")))
}

fn unsigned(param: &Param, key: &str) -> Result<usize> {
    usize::try_from(param.value(key)?.to_i64()?)
        .map_err(|_| bad(&format!("'{key}' must not be negative")))
}

impl AScore {
    /// The source defaults, in the layout of the source `defaults_`.
    ///
    /// # Errors
    ///
    /// Never in practice; as [`AScore::to_param`].
    pub fn defaults() -> Result<Param> {
        Self::default().to_param()
    }

    /// The parameters of this instance in the layout of [`AScore::defaults`],
    /// carrying the current values (source `getParameters`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `max_peptide_length` or
    /// `max_permutations` does not fit the source's 32-bit integer parameter.
    pub fn to_param(&self) -> Result<Param> {
        let mut tolerance = entry(
            "fragment_mass_tolerance",
            ParamValue::Float(self.fragment_mass_tolerance),
            "Fragment mass tolerance for spectrum comparisons",
            false,
        );
        tolerance.min_float = 0.0;
        let mut unit = entry(
            "fragment_mass_unit",
            ParamValue::String(
                if self.fragment_tolerance_ppm {
                    "ppm"
                } else {
                    "Da"
                }
                .into(),
            ),
            "Unit of fragment mass tolerance",
            false,
        );
        unit.valid_strings = vec!["Da".into(), "ppm".into()];
        let mut length = entry(
            "max_peptide_length",
            ParamValue::Integer(source_int(self.max_peptide_length, "max_peptide_length")?),
            "Restrict scoring to peptides with a length no greater than this value ('0' for 'no restriction')",
            true,
        );
        length.min_int = 0;
        let mut permutations = entry(
            "max_num_perm",
            ParamValue::Integer(source_int(self.max_permutations, "max_num_perm")?),
            "Maximum number of permutations a sequence can have to be processed ('0' for 'no restriction')",
            true,
        );
        permutations.min_int = 0;
        let unambiguous = entry(
            "unambiguous_score",
            ParamValue::Integer(i64::from(self.unambiguous_score)),
            "Score to use for unambiguous assignments, where all sites on a peptide are phosphorylated. (Note: If a peptide is not phosphorylated at all, its score is set to '-1'.)",
            true,
        );
        Param::from_root(ParamNode {
            name: String::new(),
            description: String::new(),
            entries: vec![tolerance, unit, length, permutations, unambiguous],
            nodes: Vec::new(),
        })
    }

    /// Build an AScore from parameters, as source `setParameters` followed by
    /// `updateMembers_`.
    ///
    /// Missing parameters take their defaults; unknown parameters are ignored,
    /// as the source only warns (see [`AScore::from_param_with_warnings`]).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a parameter has the wrong value type
    /// or violates its restriction (the source `Exception::InvalidParameter`),
    /// or when `unambiguous_score` does not fit the source's 32-bit integer.
    pub fn from_param(parameters: &Param) -> Result<Self> {
        Ok(Self::from_param_with_warnings(parameters)?.0)
    }

    /// As [`AScore::from_param`], also returning the unknown-parameter
    /// warnings the source logs.
    ///
    /// # Errors
    ///
    /// As [`AScore::from_param`].
    pub fn from_param_with_warnings(parameters: &Param) -> Result<(Self, Vec<String>)> {
        let mut handler = DefaultParamHandler::new(ASCORE_NAME)?;
        handler.set_defaults(Self::defaults()?)?;
        handler.set_parameters_with(parameters, Self::from_complete_param)
    }

    fn from_complete_param(param: &Param) -> Result<Self> {
        let unambiguous_score = i32::try_from(param.value("unambiguous_score")?.to_i64()?)
            .map_err(|_| bad("'unambiguous_score' does not fit the source 32-bit parameter"))?;
        Ok(Self {
            fragment_mass_tolerance: param.value("fragment_mass_tolerance")?.to_f64()?,
            fragment_tolerance_ppm: param.value("fragment_mass_unit")?.as_str()? == "ppm",
            max_peptide_length: unsigned(param, "max_peptide_length")?,
            max_permutations: unsigned(param, "max_num_perm")?,
            unambiguous_score,
        })
    }

    fn tolerance(&self) -> Result<Tolerance> {
        if !self.fragment_mass_tolerance.is_finite() || self.fragment_mass_tolerance < 0.0 {
            return Err(bad(
                "AScore fragment mass tolerance must be finite and not negative",
            ));
        }
        Ok(if self.fragment_tolerance_ppm {
            Tolerance::Ppm(self.fragment_mass_tolerance)
        } else {
            Tolerance::Absolute(self.fragment_mass_tolerance)
        })
    }

    /// Compute the AScore of a peptide hit and return the rescored hit.
    ///
    /// The returned hit is a copy of `hit` with:
    ///
    /// * the sequence of the best-scoring phosphosite placement. Phospho
    ///   information is moved; every other modification is kept as it was.
    ///   (The source's brief says other modifications are dropped; its code
    ///   keeps them, and so does this port.)
    /// * metadata `search_engine_sequence`, the input sequence (the source's
    ///   `@note` names it `Search_engine_sequence`; the key it writes, and this
    ///   port writes, is lowercase);
    /// * metadata `AScore_pep_score`, the best weighted placement score;
    /// * metadata `AScore_1`, `AScore_2`, ... one per phosphorylation of the
    ///   best placement, in residue order;
    /// * as score, the lowest of the `AScore_<n>` values.
    ///
    /// The score is reset to `-1` and nothing else changes when the spectrum is
    /// empty, the peptide has no phosphorylation or no S/T/Y, the peptide is
    /// longer than a positive `max_peptide_length`, or there are more
    /// placements than a positive `max_permutations`. When every S/T/Y is
    /// phosphorylated the score is `unambiguous_score` and nothing else changes.
    /// When the best two placements score equally, every `AScore_<n>` is `0`.
    ///
    /// `real_spectrum` is sorted by m/z in place when it is not sorted, as the
    /// source does through its in/out reference.
    ///
    /// A peptide carrying more phosphorylations than it has S/T/Y residues (a
    /// phospho on another residue) has no placement to score; the source reads
    /// out of bounds there, and this port returns the reset hit (score `-1`)
    /// as for the other aborted cases.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the tolerance is negative or not
    /// finite, when the placement count exceeds [`MAX_ASCORE_PERMUTATIONS`],
    /// when the spectrum has a negative m/z or spans more than
    /// [`MAX_ASCORE_WINDOWS`] windows, when a match probability exceeds one or
    /// a matched-ion count exceeds its ion count (both throw
    /// `std::invalid_argument` in the source's binomial), or when a score is not
    /// finite and so cannot be stored as metadata. Errors from sequence parsing,
    /// spectrum generation and sorting are passed on. On error `hit` and the
    /// peaks of `real_spectrum` are unchanged, except that the spectrum may
    /// already be sorted.
    pub fn compute(&self, hit: &PeptideHit, real_spectrum: &mut MSSpectrum) -> Result<PeptideHit> {
        self.tolerance()?;
        let mut phospho = hit.clone();
        phospho.score = -1.0;
        if real_spectrum.is_empty() {
            return Ok(phospho);
        }
        let sequence_str = phospho.sequence.to_string();
        let unmodified_sequence_str = phospho.sequence.as_str().to_owned();

        let number_of_phosphorylation_events = Self::number_of_phospho_events(&sequence_str);
        let seq_without_phospho = Self::remove_phosphosites_from_sequence(&sequence_str)?;

        if self.max_peptide_length > 0 && unmodified_sequence_str.len() > self.max_peptide_length {
            return Ok(phospho);
        }

        let sites = Self::get_sites(&unmodified_sequence_str);
        let number_of_sty = sites.len();
        if number_of_phosphorylation_events == 0 || number_of_sty == 0 {
            return Ok(phospho);
        }
        if number_of_sty == number_of_phosphorylation_events {
            phospho.score = f64::from(self.unambiguous_score);
            return Ok(phospho);
        }
        if number_of_phosphorylation_events > number_of_sty {
            return Ok(phospho);
        }
        let placements = Self::permutation_count(number_of_sty, number_of_phosphorylation_events);
        if self.max_permutations > 0 && placements > self.max_permutations {
            return Ok(phospho);
        }
        let permutations = Self::compute_permutations(&sites, number_of_phosphorylation_events)?;

        let th_spectra = self.create_theoretical_spectra(&permutations, &seq_without_phospho)?;

        if !real_spectrum.is_sorted() {
            real_spectrum.sort_by_position()?;
        }
        let windows_top10 = Self::peak_picking_per_windows_in_spectrum(real_spectrum)?;

        let last_mz = real_spectrum.peaks[real_spectrum.len() - 1].mz;
        let base_match_probability = self.compute_base_probability(last_mz);

        let peptide_site_scores = self.calculate_permutation_peptide_scores(
            &th_spectra,
            &windows_top10,
            base_match_probability,
        )?;
        let ranking = Self::rank_weighted_permutation_peptide_scores(&peptide_site_scores)?;

        // At least two placements exist: there are more S/T/Y than events.
        let (peptide1_score, best) = ranking[ranking.len() - 1];
        let (peptide2_score, _) = ranking[ranking.len() - 2];
        phospho.sequence = AASequence::parse(&th_spectra[best].name)?;
        phospho.metadata.insert(
            "search_engine_sequence".into(),
            MetaValue::from(hit.sequence.to_string()),
        );
        phospho.metadata.insert(
            "AScore_pep_score".into(),
            MetaValue::try_from(peptide1_score)?,
        );

        let phospho_sites = Self::determine_highest_scoring_permutations(
            &peptide_site_scores,
            &permutations,
            &ranking,
        )?;

        let mut best_ascore = f64::MAX; // the lower the better
        for (index, site) in phospho_sites.iter().enumerate() {
            let mut ascore = 0.0;
            if peptide1_score != peptide2_score {
                let [first, second] = self.compute_site_determining_ions(&th_spectra, site)?;
                let p = site.peak_depth as f64 * base_match_probability;
                let mut n_first = 0;
                for window in &windows_top10 {
                    n_first += self.number_of_matched_ions(&first, window, site.peak_depth)?;
                }
                let p_first = Self::compute_cumulative_score(first.len(), n_first, p)?;
                let mut n_second = 0;
                for window in &windows_top10 {
                    n_second += self.number_of_matched_ions(&second, window, site.peak_depth)?;
                }
                let p_second = Self::compute_cumulative_score(second.len(), n_second, p)?;
                // abs avoids -0 scores, as in the source.
                let score_first = (-10.0 * p_first.log10()).abs();
                let score_second = (-10.0 * p_second.log10()).abs();
                ascore = score_first - score_second;
            }
            if ascore < best_ascore {
                best_ascore = ascore;
            }
            phospho.metadata.insert(
                format!("AScore_{}", index + 1),
                MetaValue::try_from(ascore)?,
            );
        }
        phospho.score = best_ascore;
        Ok(phospho)
    }

    /// Compare two m/z values within the fragment tolerance (source
    /// `compareMZ_`).
    ///
    /// `Less` when `mz1 - mz2 < -tolerance`, `Greater` when it exceeds
    /// `tolerance`, `Equal` otherwise. In ppm mode the tolerance is taken
    /// relative to the mean of the two values. The source returns `-1`, `1`
    /// and `0`.
    pub fn compare_mz(&self, mz1: f64, mz2: f64) -> Ordering {
        let mut tolerance = self.fragment_mass_tolerance;
        let error = mz1 - mz2;
        if self.fragment_tolerance_ppm {
            let avg_mass = (mz1 + mz2) / 2.0;
            tolerance = tolerance * avg_mass / 1e6;
        }
        if error < -tolerance {
            Ordering::Less
        } else if error > tolerance {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }

    /// The peaks of `first` that have no counterpart in `second` (source
    /// `getSpectrumDifference_`).
    ///
    /// Works like `std::set_difference` with [`AScore::compare_mz`] as the
    /// comparison; the source reimplements it because a tolerance comparison
    /// is not a total order. When two peaks match, neither is output, and every
    /// further peak of `second` matching the `first` peak and every further
    /// peak of `first` matching the `second` peak is skipped too: several ions
    /// within one tolerance cannot be site determining. Both inputs are
    /// expected in ascending m/z.
    ///
    /// The source's skip loops dereference the end iterator after stepping
    /// onto it; the value read is discarded because the loop condition then
    /// fails, so skipping up to the end, as this port does, is the same result.
    pub fn spectrum_difference(&self, first: &[Peak1D], second: &[Peak1D]) -> Vec<Peak1D> {
        let mut result = Vec::new();
        let (mut i, mut j) = (0, 0);
        while i < first.len() && j < second.len() {
            let mz1 = first[i].mz;
            let mz2 = second[j].mz;
            match self.compare_mz(mz1, mz2) {
                Ordering::Less => {
                    result.push(first[i]);
                    i += 1;
                }
                Ordering::Greater => j += 1,
                Ordering::Equal => {
                    j += 1;
                    while j < second.len() && self.compare_mz(mz1, second[j].mz) == Ordering::Equal
                    {
                        j += 1;
                    }
                    i += 1;
                    while i < first.len() && self.compare_mz(first[i].mz, mz2) == Ordering::Equal {
                        i += 1;
                    }
                }
            }
        }
        result.extend_from_slice(&first[i..]);
        result
    }

    /// The site-determining ions of a site and its runner-up (source
    /// `computeSiteDeterminingIons_`).
    ///
    /// Element 0 holds the ions of placement `candidates.seq_1` that are not in
    /// placement `candidates.seq_2`, element 1 the reverse, each from
    /// [`AScore::spectrum_difference`] and sorted by m/z. The spectra carry
    /// only the peaks.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either placement index is out of
    /// range of `th_spectra` (the source reads out of bounds).
    pub fn compute_site_determining_ions(
        &self,
        th_spectra: &[MSSpectrum],
        candidates: &ProbablePhosphoSites,
    ) -> Result<[MSSpectrum; 2]> {
        let spectrum_first = th_spectra
            .get(candidates.seq_1)
            .ok_or_else(|| bad("AScore candidate seq_1 is out of range"))?;
        let spectrum_second = th_spectra
            .get(candidates.seq_2)
            .ok_or_else(|| bad("AScore candidate seq_2 is out of range"))?;
        let mut first = MSSpectrum {
            peaks: self.spectrum_difference(&spectrum_first.peaks, &spectrum_second.peaks),
            ..MSSpectrum::default()
        };
        let mut second = MSSpectrum {
            peaks: self.spectrum_difference(&spectrum_second.peaks, &spectrum_first.peaks),
            ..MSSpectrum::default()
        };
        first.sort_by_position()?;
        second.sort_by_position()?;
        Ok([first, second])
    }

    /// All possible phosphorylation sites (source `getSites_`): the byte
    /// positions of `S`, `T` and `Y` in an unmodified sequence, ascending.
    pub fn get_sites(unmodified_sequence: &str) -> Vec<usize> {
        unmodified_sequence
            .bytes()
            .enumerate()
            .filter(|&(_, residue)| matches!(residue, b'S' | b'T' | b'Y'))
            .map(|(index, _)| index)
            .collect()
    }

    fn permutation_count(sites: usize, events: usize) -> usize {
        if events > sites {
            return 0;
        }
        let events = events.min(sites - events);
        let mut value: usize = 1;
        for step in 0..events {
            // value * (sites - step) / (step + 1) is always integral; saturate
            // rather than overflow, since anything that large is refused.
            value = match value.checked_mul(sites - step) {
                Some(product) => product / (step + 1),
                None => return usize::MAX,
            };
        }
        value
    }

    /// Every set of `n_phosphorylation_events` sites (source
    /// `computePermutations_`).
    ///
    /// The sets keep the order of `sites` and are listed in lexicographic order
    /// of their positions in `sites`, as the source's recursion (all sets with
    /// the first site, then all without it) produces. Zero events give no set;
    /// as many events as sites give the one set `sites` itself.
    ///
    /// More events than sites give no set. The source recurses into
    /// `sites[0]` of an empty vector there; the port returns the empty
    /// enumeration, which is the mathematical answer. The source's `Int`
    /// argument cannot be negative here.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there would be more than
    /// [`MAX_ASCORE_PERMUTATIONS`] sets; the count is checked before anything
    /// is allocated.
    pub fn compute_permutations(
        sites: &[usize],
        n_phosphorylation_events: usize,
    ) -> Result<Vec<Vec<usize>>> {
        let k = n_phosphorylation_events;
        if k == 0 || k > sites.len() {
            return Ok(Vec::new());
        }
        let total = Self::permutation_count(sites.len(), k);
        if total > MAX_ASCORE_PERMUTATIONS {
            return Err(bad(
                "AScore placement count exceeds MAX_ASCORE_PERMUTATIONS",
            ));
        }
        let mut permutations = Vec::with_capacity(total);
        let mut indices: Vec<usize> = (0..k).collect();
        loop {
            permutations.push(indices.iter().map(|&index| sites[index]).collect());
            // Advance to the next combination in lexicographic order.
            let mut position = k;
            let mut advanced = false;
            while position > 0 {
                position -= 1;
                if indices[position] != position + sites.len() - k {
                    indices[position] += 1;
                    for next in position + 1..k {
                        indices[next] = indices[next - 1] + 1;
                    }
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
        }
        Ok(permutations)
    }

    /// Matched ions between a theoretical spectrum and the top `depth` peaks of
    /// one window (source `numberOfMatchedIons_`).
    ///
    /// `window` is expected in descending intensity, as
    /// [`AScore::peak_picking_per_windows_in_spectrum`] returns it: its first
    /// `depth` peaks are kept, sorted by m/z, and the theoretical peaks paired
    /// with one of them by the source's `MatchedIterator` walk are counted
    /// (one per theoretical peak, inclusive tolerance in single precision, ppm
    /// taken from the theoretical m/z).
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsortedData`] when `th` is not sorted by m/z (the
    /// source requires it and does not check), and [`Error::InvalidValue`] when
    /// the tolerance is negative, not finite or outside the `f32` range.
    pub fn number_of_matched_ions(
        &self,
        th: &MSSpectrum,
        window: &MSSpectrum,
        depth: usize,
    ) -> Result<usize> {
        let tolerance = self.tolerance()?;
        if !th.is_sorted() {
            return Err(Error::UnsortedData);
        }
        let mut window_reduced = MSSpectrum {
            peaks: window.peaks.iter().take(depth).copied().collect(),
            ..MSSpectrum::default()
        };
        window_reduced.sort_by_position()?;
        Ok(matched_alignment(th, &window_reduced, tolerance)?.len())
    }

    /// The weighted peptide score of Beausoleil et al., page 1291 (source
    /// `peptideScore_`).
    ///
    /// `(0.5 s1 + 0.75 s2 + s3 + s4 + s5 + s6 + 0.75 s7 + 0.5 s8 + 0.25 s9 +
    /// 0.25 s10) / 7` over the scores at peak depths 1 to 10, summed in that
    /// order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] unless there are exactly
    /// [`ASCORE_PEAK_DEPTHS`] scores; the source states this as a debug-only
    /// precondition and reads out of bounds without it.
    pub fn peptide_score(scores: &[f64]) -> Result<f64> {
        if scores.len() != ASCORE_PEAK_DEPTHS {
            return Err(bad(
                "AScore peptide score needs a score for every peak level",
            ));
        }
        Ok((scores[0] * 0.5
            + scores[1] * 0.75
            + scores[2]
            + scores[3]
            + scores[4]
            + scores[5]
            + scores[6] * 0.75
            + scores[7] * 0.5
            + scores[8] * 0.25
            + scores[9] * 0.25)
            / 7.0)
    }

    /// The runner-up placement and best peak depth for every site of the best
    /// placement (source `determineHighestScoringPermutations_`).
    ///
    /// `ranking` is the ascending weighted ranking from
    /// [`AScore::rank_weighted_permutation_peptide_scores`]; its last entry is
    /// the best placement. For each site of the best placement, in placement
    /// order:
    ///
    /// 1. the next placement down the ranking that lacks this site but has
    ///    every other site of the best placement becomes `seq_2`, and the site
    ///    it uses instead becomes `second`;
    /// 2. the peak depth (1 to 10) at which the unweighted score of the best
    ///    placement exceeds that of `seq_2` by the most becomes `peak_depth`;
    ///    ties keep the lower depth and no positive difference keeps depth 1.
    ///
    /// The source sizes the result by `permutations[0]` and assumes all
    /// placements have the same size, which [`AScore::compute_permutations`]
    /// guarantees. The source notes that it assumes more placements than
    /// phosphorylations.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `ranking` or `permutations` is
    /// empty, when an index is out of range, when no lower-ranked placement
    /// qualifies as runner-up (the source walks past the start of its map),
    /// or when a score row is shorter than the best placement's row.
    pub fn determine_highest_scoring_permutations(
        peptide_site_scores: &[Vec<f64>],
        permutations: &[Vec<usize>],
        ranking: &[(f64, usize)],
    ) -> Result<Vec<ProbablePhosphoSites>> {
        let out_of_range = || bad("AScore ranking index is out of range");
        let first_permutation = permutations
            .first()
            .ok_or_else(|| bad("AScore needs at least one placement"))?;
        let &(_, best) = ranking
            .last()
            .ok_or_else(|| bad("AScore ranking must not be empty"))?;
        let best_peptide_sites = permutations.get(best).ok_or_else(out_of_range)?;
        let mut sites = vec![ProbablePhosphoSites::default(); first_permutation.len()];
        if sites.len() < best_peptide_sites.len() {
            return Err(bad("AScore placements must all have the same size"));
        }

        for (i, &site) in best_peptide_sites.iter().enumerate() {
            sites[i].first = site;
            sites[i].seq_1 = best;
            let mut runner_up = None;
            for &(_, candidate) in ranking[..ranking.len() - 1].iter().rev() {
                let placement = permutations.get(candidate).ok_or_else(out_of_range)?;
                let qualifies = best_peptide_sites.iter().enumerate().all(|(j, other)| {
                    let present = placement.contains(other);
                    if j == i { !present } else { present }
                });
                if qualifies {
                    runner_up = Some(candidate);
                    break;
                }
            }
            let seq_2 =
                runner_up.ok_or_else(|| bad("AScore found no runner-up placement for a site"))?;
            sites[i].seq_2 = seq_2;
            // The site seq_2 uses that the best placement does not.
            if let Some(&other) = permutations[seq_2]
                .iter()
                .find(|position| !best_peptide_sites.contains(position))
            {
                sites[i].second = other;
            }
        }

        for site in &mut sites {
            let first_scores = peptide_site_scores
                .get(site.seq_1)
                .ok_or_else(out_of_range)?;
            let second_scores = peptide_site_scores
                .get(site.seq_2)
                .ok_or_else(out_of_range)?;
            if first_scores.len() < second_scores.len() {
                return Err(bad("AScore score rows must have the same length"));
            }
            let mut maximum_score_difference = 0.0;
            site.peak_depth = 1;
            for (depth, (phospho_at_site, no_phospho_at_site)) in
                first_scores.iter().zip(second_scores).enumerate()
            {
                let score_difference = phospho_at_site - no_phospho_at_site;
                if score_difference > maximum_score_difference {
                    maximum_score_difference = score_difference;
                    site.peak_depth = depth + 1;
                }
            }
        }
        Ok(sites)
    }

    /// The probability of a random match at a peak depth of one (source
    /// `computeBaseProbability_`).
    ///
    /// `2 * tolerance / 100`: the tolerance window's share of a 100 Th window.
    /// In ppm mode the window is converted to Th at `ppm_reference_mz`, which
    /// [`AScore::compute`] takes as the spectrum's highest m/z.
    pub fn compute_base_probability(&self, ppm_reference_mz: f64) -> f64 {
        let mut base_match_probability = 2.0 * self.fragment_mass_tolerance / 100.0;
        if self.fragment_tolerance_ppm {
            base_match_probability *= ppm_reference_mz * 1e-6;
        }
        base_match_probability
    }

    /// The cumulative binomial probability of at least `n_matched` matches in
    /// `n_trials` trials at probability `p` (source `computeCumulativeScore_`,
    /// whose arguments are `N`, `n` and `p`).
    ///
    /// Delegates to [`binomial_cdf_complement`], as the source delegates to
    /// `Math::binomial_cdf_complement`. Shared with [`PScore`](super::pscore::PScore),
    /// which the source makes a friend for exactly this function.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `n_matched > n_trials` or `p` is
    /// outside `[0, 1]`
    /// (debug-only preconditions in the source, whose binomial then throws
    /// `std::invalid_argument`), or when a count exceeds `u32::MAX`.
    pub fn compute_cumulative_score(n_trials: usize, n_matched: usize, p: f64) -> Result<f64> {
        binomial_cdf_complement(count(n_trials)?, count(n_matched)?, p)
    }

    /// The number of phosphorylations in a sequence string (source
    /// `numberOfPhosphoEvents_`): the non-overlapping occurrences of
    /// `(Phospho)`.
    pub fn number_of_phospho_events(sequence: &str) -> usize {
        sequence.matches("(Phospho)").count()
    }

    /// The peptide with every phosphorylation removed (source
    /// `removePhosphositesFromSequence_`): `(Phospho)` is deleted from the
    /// sequence string and the rest parsed again.
    ///
    /// # Errors
    ///
    /// Errors from [`AASequence::parse`] are passed on.
    pub fn remove_phosphosites_from_sequence(sequence: &str) -> Result<AASequence> {
        AASequence::parse(&sequence.replace("(Phospho)", ""))
    }

    /// One theoretical spectrum per placement (source
    /// `createTheoreticalSpectra_`).
    ///
    /// Each placement's residues get `Phospho` on `seq_without_phospho`, and a
    /// singly charged b/y spectrum is generated with the default
    /// [`TheoreticalSpectrumGenerator`]. The spectrum's name is the modified
    /// sequence's string form, which [`AScore::compute`] parses back.
    ///
    /// A placement is consumed in order while walking the residues, as in the
    /// source: a position that is not ascending, or beyond the sequence, ends
    /// the placement silently. An empty placement adds no modification (the
    /// source reads out of bounds there).
    ///
    /// # Errors
    ///
    /// Errors from [`AASequence::set_modification`] (for example `Phospho` on
    /// a residue it cannot modify) and from spectrum generation are passed on.
    pub fn create_theoretical_spectra(
        &self,
        permutations: &[Vec<usize>],
        seq_without_phospho: &AASequence,
    ) -> Result<Vec<MSSpectrum>> {
        let generator = TheoreticalSpectrumGenerator::default();
        let mut th_spectra = Vec::with_capacity(permutations.len());
        for permutation in permutations {
            let mut seq = seq_without_phospho.clone();
            let mut permu = 0;
            for residue in 0..seq.len() {
                if permu == permutation.len() {
                    break;
                }
                if residue == permutation[permu] {
                    seq.set_modification(residue, "Phospho")?;
                    permu += 1;
                }
            }
            let mut spectrum = generator.generate(&seq, 1, 1, None)?;
            spectrum.name = seq.to_string();
            th_spectra.push(spectrum);
        }
        Ok(th_spectra)
    }

    /// The top 10 peaks of every 100 Th window (source
    /// `peakPickingPerWindowsInSpectrum_`).
    ///
    /// The windows run from `floor(first m/z / 100) * 100` to
    /// `ceil(last m/z / 100) * 100`; a window holds the peaks up to and
    /// including its upper bound. Each window's peaks come in descending
    /// intensity (stable), at most 10, as spectra carrying only the peaks.
    /// `real_spectrum` must be sorted by m/z, which [`AScore::compute`] ensures.
    /// An empty spectrum gives no window (the source reads its first and last
    /// peak unchecked).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when an m/z is negative (the source
    /// converts the window bound to an unsigned integer) or when the spectrum
    /// spans more than [`MAX_ASCORE_WINDOWS`] windows.
    pub fn peak_picking_per_windows_in_spectrum(
        real_spectrum: &MSSpectrum,
    ) -> Result<Vec<MSSpectrum>> {
        let (Some(front), Some(back)) = (real_spectrum.peaks.first(), real_spectrum.peaks.last())
        else {
            return Ok(Vec::new());
        };
        if !(front.mz >= 0.0 && back.mz.is_finite()) {
            return Err(bad("AScore windows need finite, non-negative m/z"));
        }
        let lower_bound = (front.mz / 100.0).floor() * 100.0;
        let upper_bound = (back.mz / 100.0).ceil() * 100.0;
        let windows = ((upper_bound - lower_bound) / 100.0).ceil();
        if !(windows >= 0.0 && windows <= MAX_ASCORE_WINDOWS as f64) {
            return Err(bad(
                "AScore spectrum spans more than MAX_ASCORE_WINDOWS windows",
            ));
        }
        let number_of_windows = windows as usize;
        let mut windows_top10 = Vec::with_capacity(number_of_windows);
        let mut current = 0;
        let mut window_upper_bound = lower_bound + 100.0;
        for _ in 0..number_of_windows {
            let mut real_window = MSSpectrum::default();
            while current < real_spectrum.len()
                && real_spectrum.peaks[current].mz <= window_upper_bound
            {
                real_window.peaks.push(real_spectrum.peaks[current]);
                current += 1;
            }
            real_window.sort_by_intensity(true)?;
            real_window.peaks.truncate(10);
            windows_top10.push(real_window);
            window_upper_bound += 100.0;
        }
        Ok(windows_top10)
    }

    /// The unweighted scores of every placement at peak depths 1 to 10
    /// (source `calculatePermutationPeptideScores_`, Beausoleil et al.
    /// Figure 3b).
    ///
    /// For depth `i`, the ions of the placement's theoretical spectrum matched
    /// by the top `i` peaks of every window are summed into `n`, and the score
    /// is `|-10 log10 P(X >= n)|` for `X ~ B(N, i * base_match_probability)`
    /// with `N` the theoretical ion count. The absolute value avoids `-0`, as
    /// in the source. `base_match_probability` is the source's
    /// `base_match_probability_` member, from [`AScore::compute_base_probability`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a match probability exceeds one or
    /// `n` exceeds `N` (an ion matching in two adjacent windows is counted in
    /// both); the source's binomial throws `std::invalid_argument` for both.
    /// Errors from [`AScore::number_of_matched_ions`] are passed on.
    pub fn calculate_permutation_peptide_scores(
        &self,
        th_spectra: &[MSSpectrum],
        windows_top10: &[MSSpectrum],
        base_match_probability: f64,
    ) -> Result<Vec<Vec<f64>>> {
        let mut permutation_peptide_scores = Vec::with_capacity(th_spectra.len());
        for th in th_spectra {
            let n_trials = th.len();
            let mut site_score = Vec::with_capacity(ASCORE_PEAK_DEPTHS);
            for depth in 1..=ASCORE_PEAK_DEPTHS {
                let mut n = 0;
                for window in windows_top10 {
                    n += self.number_of_matched_ions(th, window, depth)?;
                }
                let p = depth as f64 * base_match_probability;
                let cumulative_score = Self::compute_cumulative_score(n_trials, n, p)?;
                site_score.push((-10.0 * cumulative_score.log10()).abs());
            }
            permutation_peptide_scores.push(site_score);
        }
        Ok(permutation_peptide_scores)
    }

    /// The placements ranked by weighted score, ascending (source
    /// `rankWeightedPermutationPeptideScores_`).
    ///
    /// Each entry is `(peptide_score, placement index)`, with
    /// [`AScore::peptide_score`] of the placement's row. The source returns a
    /// `std::multimap<double, Size>`, in which equal scores keep insertion
    /// order; this sorted vector keeps it too (a stable sort), so the last
    /// entry is the best placement and, among equal best scores, the one with
    /// the highest index, as the source's `rbegin()` gives.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a row does not have
    /// [`ASCORE_PEAK_DEPTHS`] scores, or when a weighted score is NaN (the
    /// source's map would be ordered inconsistently).
    pub fn rank_weighted_permutation_peptide_scores(
        peptide_site_scores: &[Vec<f64>],
    ) -> Result<Vec<(f64, usize)>> {
        let mut ranking = Vec::with_capacity(peptide_site_scores.len());
        for (index, scores) in peptide_site_scores.iter().enumerate() {
            let weighted_score = Self::peptide_score(scores)?;
            if weighted_score.is_nan() {
                return Err(bad("AScore weighted score must not be NaN"));
            }
            ranking.push((weighted_score, index));
        }
        // NaN is refused above, so partial_cmp is total here; it also keeps
        // -0 and +0 equal, as the source's `<` does.
        ranking.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
        Ok(ranking)
    }
}
