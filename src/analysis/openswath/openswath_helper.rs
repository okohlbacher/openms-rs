// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Helpers used by several OpenSWATH tools: precursor ids, SWATH-window
//! transition selection, SWATH-map sanity checks, the library RT range and the
//! best feature per transition group.
//!
//! Ports `ANALYSIS/OPENSWATH/OpenSwathHelper.h` with `OpenSwathHelper.cpp` (a
//! class of static members, which become free functions). See
//! [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/OPENSWATH_DATAACCESS_SUPPORT.md).
//!
//! Both `selectSwathTransitions` overloads (the `OpenMS::TargetedExperiment`
//! one is [`select_swath_transitions_targeted`](crate::analysis::openswath::openswath_helper::select_swath_transitions_targeted)) and the template
//! `checkSwathMapAndSelectTransitions`, whose two instantiations are one generic
//! function over [`SwathSelectable`](crate::analysis::openswath::openswath_helper::SwathSelectable), are ported, and so is `sampleExperiment`
//! ([`sample_experiment`](crate::analysis::openswath::openswath_helper::sample_experiment)), which draws through
//! [`RandomShuffler`](crate::math::random_shuffler::RandomShuffler).

use super::swath_map::SwathMap;
use super::targeted_experiment::TargetedExperiment;
use super::transition_experiment::{LightCompound, LightTargetedExperiment, LightTransition};
use crate::kernel::MSExperiment;
use crate::kernel::mrm::MRMTransitionGroup;
use crate::math::random_shuffler::RandomShuffler;
use crate::math::source_sort::source_sort_by;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// The unique precursor id `"<transition_group_id>_Precursor_i<isotope>"`; the
/// monoisotopic precursor of `"groupID"` is `"groupID_Precursor_i0"`.
pub fn compute_precursor_id(transition_group_id: &str, isotope: i32) -> String {
    format!("{transition_group_id}_Precursor_i{isotope}")
}

/// The transition group (peptide or compound) id of a precursor id, reversing
/// [`compute_precursor_id`].
///
/// The id is split at every `_`: with exactly three parts the first is
/// returned; with more, all but the last two are joined with `_` again (so a
/// group id may itself contain underscores, even doubled ones); with fewer, or
/// for the empty string, the result is empty.
pub fn compute_transition_group_id(precursor_id: &str) -> String {
    if precursor_id.is_empty() {
        return String::new();
    }
    let parts: Vec<&str> = precursor_id.split('_').collect();
    match parts.len() {
        3 => parts[0].to_string(),
        n if n > 3 => parts[..n - 2].join("_"),
        _ => String::new(),
    }
}

/// Whether a precursor m/z lies inside a SWATH window by the source rule:
/// strictly between `lower` and `upper`, and at least `min_upper_edge_dist`
/// below the upper edge.
fn in_window(precursor_mz: f64, lower: f64, upper: f64, min_upper_edge_dist: f64) -> bool {
    lower < precursor_mz
        && precursor_mz < upper
        && (upper - precursor_mz).abs() >= min_upper_edge_dist
}

/// Append to `selected` the transitions of `targeted_exp` whose precursor lies
/// in the SWATH window `(lower, upper)` (both exclusive) and at least
/// `min_upper_edge_dist` Th below the upper edge, then the compounds those
/// transitions reference and the proteins those compounds reference, each in
/// library order.
///
/// `selected` is appended to, not cleared, as in the source.
pub fn select_swath_transitions(
    targeted_exp: &LightTargetedExperiment,
    selected: &mut LightTargetedExperiment,
    min_upper_edge_dist: f64,
    lower: f64,
    upper: f64,
) {
    let mut matching_compounds = BTreeSet::new();
    for tr in &targeted_exp.transitions {
        if in_window(tr.precursor_mz, lower, upper, min_upper_edge_dist) {
            selected.transitions.push(tr.clone());
            matching_compounds.insert(tr.peptide_ref.as_str());
        }
    }
    let mut matching_proteins = BTreeSet::new();
    for compound in &targeted_exp.compounds {
        if matching_compounds.contains(compound.id.as_str()) {
            selected.compounds.push(compound.clone());
            matching_proteins.extend(compound.protein_refs.iter().map(String::as_str));
        }
    }
    for protein in &targeted_exp.proteins {
        if matching_proteins.contains(protein.id.as_str()) {
            selected.proteins.push(protein.clone());
        }
    }
}

/// The `OpenMS::TargetedExperiment` overload of [`select_swath_transitions`]
/// (source `selectSwathTransitions(const TargetedExperiment&,
/// TargetedExperiment&, ...)`).
///
/// Unlike the light overload this one *replaces* the peptides and proteins of
/// `selected` by all of those of `targeted_exp` (source `setPeptides`,
/// `setProteins`) and only then appends the transitions whose precursor lies
/// in the window by the rule of [`select_swath_transitions`]; compounds are not
/// filtered.
pub fn select_swath_transitions_targeted(
    targeted_exp: &TargetedExperiment,
    selected: &mut TargetedExperiment,
    min_upper_edge_dist: f64,
    lower: f64,
    upper: f64,
) {
    selected.set_peptides(targeted_exp.peptides().to_vec());
    selected.set_proteins(targeted_exp.proteins().to_vec());
    for tr in targeted_exp.transitions() {
        if in_window(tr.precursor_mz(), lower, upper, min_upper_edge_dist) {
            selected.add_transition(tr.clone());
        }
    }
}

/// The two transition-list types `checkSwathMapAndSelectTransitions` is
/// instantiated for: the light experiment and the heavyweight
/// `TargetedExperiment`.
pub trait SwathSelectable: Sized {
    /// `selectSwathTransitions` for this type.
    fn select_swath_transitions(
        &self,
        selected: &mut Self,
        min_upper_edge_dist: f64,
        lower: f64,
        upper: f64,
    );

    /// Whether the experiment holds no transition (`getTransitions().empty()`).
    fn has_no_transitions(&self) -> bool;
}

impl SwathSelectable for LightTargetedExperiment {
    fn select_swath_transitions(
        &self,
        selected: &mut Self,
        min_upper_edge_dist: f64,
        lower: f64,
        upper: f64,
    ) {
        select_swath_transitions(self, selected, min_upper_edge_dist, lower, upper);
    }

    fn has_no_transitions(&self) -> bool {
        self.transitions.is_empty()
    }
}

impl SwathSelectable for TargetedExperiment {
    fn select_swath_transitions(
        &self,
        selected: &mut Self,
        min_upper_edge_dist: f64,
        lower: f64,
        upper: f64,
    ) {
        select_swath_transitions_targeted(self, selected, min_upper_edge_dist, lower, upper);
    }

    fn has_no_transitions(&self) -> bool {
        self.transitions().is_empty()
    }
}

/// Match every transition with its best window across m/z and ion mobility.
///
/// For diaPASEF, several windows can share an m/z range and differ in ion
/// mobility. A window matches a transition when the precursor ion mobility lies
/// strictly inside `(im_lower, im_upper)` and the precursor m/z lies in the
/// window by the rule of [`select_swath_transitions`]. Of several matching
/// windows the one whose ion-mobility centre is strictly closer to the
/// precursor ion mobility wins; on a tie the earlier window is kept.
///
/// Returns, per transition, the index into `swath_maps` of its window, or
/// `None` (the source's -1). The source resizes the caller's map, keeping any
/// entries it already held; the port returns a fresh one.
///
/// The source states as a (debug-only) precondition that all transitions carry
/// an ion mobility, while checking only that one does; neither is checked
/// here. A transition without one (-1) matches only windows whose lower bound
/// is below -1.
pub fn select_swath_transitions_pasef(
    transition_exp: &LightTargetedExperiment,
    min_upper_edge_dist: f64,
    swath_maps: &[SwathMap],
) -> Vec<Option<usize>> {
    let mut tr_win_map: Vec<Option<usize>> = vec![None; transition_exp.transitions.len()];
    let im_center = |m: &SwathMap| (m.im_lower + m.im_upper) / 2.0;
    for (i, map) in swath_maps.iter().enumerate() {
        for (k, tr) in transition_exp.transitions.iter().enumerate() {
            let im = tr.precursor_im;
            if map.im_lower < im
                && im < map.im_upper
                && in_window(tr.precursor_mz, map.lower, map.upper, min_upper_edge_dist)
            {
                match tr_win_map[k] {
                    None => tr_win_map[k] = Some(i),
                    Some(old) => {
                        let im_old = (im_center(&swath_maps[old]) - im).abs();
                        let im_new = (im_center(map) - im).abs();
                        if im_old > im_new {
                            tr_win_map[k] = Some(i);
                        }
                    }
                }
            }
        }
    }
    tr_win_map
}

/// The lower edge, upper edge and centre (all in Th) of a SWATH map, after
/// sanity checks over the whole map.
///
/// The window is that of the first spectrum's precursor: centre `mz`, lower
/// edge `mz - lower offset`, upper edge `mz + upper offset`. Every spectrum
/// must have exactly one precursor, the MS level of the first spectrum, and a
/// precursor whose m/z and both isolation offsets are within 0.1 Th of the
/// first one's (otherwise extracting a chromatogram across them makes no sense).
///
/// # Errors
///
/// [`Error::InvalidValue`] (source `Exception::IllegalArgument`) when the map
/// has no spectra, the first spectrum has no precursor, or any check fails.
pub fn check_swath_map(swath_map: &MSExperiment) -> Result<(f64, f64, f64)> {
    let first = swath_map
        .spectra
        .first()
        .and_then(|s| s.precursors.first())
        .ok_or_else(|| Error::InvalidValue("Swath map has no Spectra".into()))?;
    let lower = first.mz - first.isolation_window_lower_offset;
    let upper = first.mz + first.isolation_window_upper_offset;
    let center = first.mz;
    let expected_ms_level = swath_map.spectra[0].ms_level;
    for (index, spectrum) in swath_map.spectra.iter().enumerate() {
        if spectrum.precursors.len() != 1 {
            return Err(Error::InvalidValue(format!(
                "Scan {index} does not have exactly one precursor."
            )));
        }
        if spectrum.ms_level != expected_ms_level {
            return Err(Error::InvalidValue(format!(
                "Scan {index} is of a different MS level than the first scan."
            )));
        }
        let prec = &spectrum.precursors[0];
        if (prec.mz - first.mz).abs() > 0.1
            || (prec.isolation_window_lower_offset - first.isolation_window_lower_offset).abs()
                > 0.1
            || (prec.isolation_window_upper_offset - first.isolation_window_upper_offset).abs()
                > 0.1
        {
            return Err(Error::InvalidValue(format!(
                "Scan {index} has a different precursor isolation window than the first scan."
            )));
        }
    }
    Ok((lower, upper, center))
}

/// Check a SWATH map and select the transitions of its window in one call.
///
/// Returns `Ok(false)` when the map has no spectra or its first spectrum has
/// no precursor (not a SWATH map), and when no transition falls into the
/// window; the source prints a warning to standard error in both cases, which
/// the caller can do on `false`. Otherwise the window is computed by
/// [`check_swath_map`], transitions are appended to `selected` by
/// [`select_swath_transitions`], and the result is `Ok(true)`.
///
/// The source checks `selected` for emptiness after appending, so a `selected`
/// that already held transitions never reports an empty window; the port keeps
/// that rule.
///
/// # Errors
///
/// Those of [`check_swath_map`].
pub fn check_swath_map_and_select_transitions<T: SwathSelectable>(
    exp: &MSExperiment,
    targeted_exp: &T,
    selected: &mut T,
    min_upper_edge_dist: f64,
) -> Result<bool> {
    if exp
        .spectra
        .first()
        .is_none_or(|spectrum| spectrum.precursors.is_empty())
    {
        return Ok(false);
    }
    let (lower, upper, _center) = check_swath_map(exp)?;
    targeted_exp.select_swath_transitions(selected, min_upper_edge_dist, lower, upper);
    Ok(!selected.has_no_transitions())
}

/// The `(min, max)` retention time over the compounds of `exp`.
///
/// Comparisons start from the first compound; a NaN RT never compares smaller
/// or larger, so it is ignored unless the first compound carries it, in which
/// case it stays (both ends NaN if no other compound compares against it), as
/// in the source.
///
/// # Errors
///
/// [`Error::InvalidValue`] when `exp` has no compounds (source
/// `Exception::IllegalArgument`, "Input list of targets is empty.").
pub fn estimate_rt_range(exp: &LightTargetedExperiment) -> Result<(f64, f64)> {
    let first = exp
        .compounds
        .first()
        .ok_or_else(|| Error::InvalidValue("Input list of targets is empty.".into()))?;
    let (mut min, mut max) = (first.rt, first.rt);
    for compound in &exp.compounds {
        if compound.rt < min {
            min = compound.rt;
        }
        if compound.rt > max {
            max = compound.rt;
        }
    }
    Ok((min, max))
}

/// Upper bound on `bins` times the number of candidate compounds, the work of
/// the two sampling passes; checked before any bucket is built.
pub const MAX_SAMPLING_WORK: usize = 1 << 28;

/// Sample a retention-time-stratified subset of a library, for iRT calibration
/// (source `sampleExperiment`).
///
/// 1. Candidates are the compounds referenced by at least one non-decoy
///    transition; those whose sequence is in `priority_peptides` form the
///    priority list, the others the regular list.
/// 2. With `sort_by_intensity`, both lists are sorted by the summed library
///    intensity of their non-decoy transitions, descending, with the crate's
///    emulation of the Release build's `std::sort` (ties are not stable), and
///    the regular list is cut to its first
///    `max(1, floor(len * top_fraction))` entries.
/// 3. Fewer than three candidates in total is an error.
/// 4. The retention-time range of all candidates is cut into `bins` equal
///    bins. Priority peptides are sampled first: per bin, the not yet picked
///    (by *sequence*) priority peptides with `rt` in `[lo, hi)` are shuffled
///    with a [`RandomShuffler`] seeded with `seed` and the first
///    `peptides_per_bin` are taken. Then, when fewer than
///    `bins * peptides_per_bin` are picked and regular candidates exist, every
///    bin is filled up to `peptides_per_bin` counting everything already picked
///    in it. The last bin's upper edge is `rt_max` and is exclusive, so a
///    compound exactly at the maximum is never picked, as in the source.
/// 5. The result holds the picked compounds in pick order, their non-decoy
///    transitions in library order and their proteins in library order.
///
/// `seed` 0 takes a random seed (the source's `std::random_device`), which the
/// port draws from the standard library's randomly keyed hasher.
///
/// # Errors
///
/// [`Error::InvalidValue`] for `bins` or `peptides_per_bin` of 0 and for a
/// `top_fraction` outside `(0, 1]` with `sort_by_intensity` (debug-only
/// preconditions in the source, which divides by zero or misbehaves without
/// them), for fewer than three candidates (the source's
/// `Exception::IllegalArgument`, same message), and when `bins` times the
/// number of candidates exceeds [`MAX_SAMPLING_WORK`].
pub fn sample_experiment(
    exp: &LightTargetedExperiment,
    bins: usize,
    peptides_per_bin: usize,
    seed: u32,
    sort_by_intensity: bool,
    top_fraction: f64,
    priority_peptides: &HashSet<String>,
) -> Result<LightTargetedExperiment> {
    if bins < 1 {
        return Err(Error::InvalidValue("bins must be >= 1".into()));
    }
    if peptides_per_bin < 1 {
        return Err(Error::InvalidValue("peptides_per_bin must be >= 1".into()));
    }
    if sort_by_intensity && !(top_fraction > 0.0 && top_fraction <= 1.0) {
        return Err(Error::InvalidValue(
            "top_fraction must be in (0,1] when sort_by_intensity is true".into(),
        ));
    }

    // 0) candidate selection: exclude decoys
    let good_ids: HashSet<&str> = exp
        .transitions
        .iter()
        .filter(|tr| !tr.flags.decoy())
        .map(|tr| tr.peptide_ref.as_str())
        .collect();
    let mut candidates: Vec<LightCompound> = Vec::new();
    let mut priority_candidates: Vec<LightCompound> = Vec::new();
    for compound in &exp.compounds {
        if good_ids.contains(compound.id.as_str()) {
            if priority_peptides.contains(&compound.sequence) {
                priority_candidates.push(compound.clone());
            } else {
                candidates.push(compound.clone());
            }
        }
    }

    // 1) optionally sort by library intensity and trim to the top fraction
    if sort_by_intensity && top_fraction > 0.0 && top_fraction <= 1.0 {
        let mut intensity_sum: HashMap<&str, f64> = HashMap::new();
        for tr in &exp.transitions {
            if !tr.flags.decoy() {
                *intensity_sum.entry(tr.peptide_ref.as_str()).or_insert(0.0) +=
                    tr.library_intensity;
            }
        }
        let sum_of = |id: &str| intensity_sum.get(id).copied().unwrap_or(0.0);
        source_sort_by(&mut candidates, |a, b| sum_of(&a.id) > sum_of(&b.id))?;
        // `static_cast<Size>(candidates.size() * top_fraction)`
        let max_keep = std::cmp::max(1, (candidates.len() as f64 * top_fraction) as usize);
        candidates.truncate(std::cmp::min(max_keep, candidates.len()));
        source_sort_by(&mut priority_candidates, |a, b| {
            sum_of(&a.id) > sum_of(&b.id)
        })?;
    }

    let total = priority_candidates.len() + candidates.len();
    if total < 3 {
        return Err(Error::InvalidValue(format!(
            "Insufficient candidates for sampling: {total} found, minimum 3 required for meaningful iRT calibration."
        )));
    }
    if bins.saturating_mul(total) > MAX_SAMPLING_WORK {
        return Err(Error::InvalidValue(format!(
            "sampling {bins} bins over {total} candidates exceeds the work limit of {MAX_SAMPLING_WORK}"
        )));
    }

    // 2) retention-time range over all candidates
    let mut rt_min = f64::MAX;
    let mut rt_max = f64::MIN;
    for compound in priority_candidates.iter().chain(candidates.iter()) {
        rt_min = rt_min.min(compound.rt);
        rt_max = rt_max.max(compound.rt);
    }
    let bin_width = (rt_max - rt_min) / bins as f64;
    let edges = |b: usize| {
        let lo = rt_min + b as f64 * bin_width;
        let hi = if b + 1 == bins {
            rt_max
        } else {
            lo + bin_width
        };
        (lo, hi)
    };

    // 3) priority peptides first, then the regular quota
    let mut picked: Vec<LightCompound> = Vec::new();
    let mut picked_sequences: HashSet<String> = HashSet::new();
    let mut shuffler = RandomShuffler::default();
    shuffler.seed(if seed == 0 {
        u64::from(random_seed())
    } else {
        u64::from(seed)
    });

    if !priority_candidates.is_empty() {
        for b in 0..bins {
            let (lo, hi) = edges(b);
            let mut bucket: Vec<LightCompound> = priority_candidates
                .iter()
                .filter(|c| c.rt >= lo && c.rt < hi && !picked_sequences.contains(&c.sequence))
                .cloned()
                .collect();
            if !bucket.is_empty() {
                shuffler.portable_random_shuffle(&mut bucket);
                let take = peptides_per_bin.min(bucket.len());
                for compound in bucket.into_iter().take(take) {
                    picked_sequences.insert(compound.sequence.clone());
                    picked.push(compound);
                }
            }
        }
    }

    let total_quota = bins.saturating_mul(peptides_per_bin);
    if picked.len() < total_quota && !candidates.is_empty() {
        for b in 0..bins {
            let (lo, hi) = edges(b);
            let priority_in_bin = picked.iter().filter(|p| p.rt >= lo && p.rt < hi).count();
            let needed = peptides_per_bin.saturating_sub(priority_in_bin);
            if needed > 0 {
                let mut bucket: Vec<LightCompound> = candidates
                    .iter()
                    .filter(|c| c.rt >= lo && c.rt < hi && !picked_sequences.contains(&c.sequence))
                    .cloned()
                    .collect();
                if !bucket.is_empty() {
                    shuffler.portable_random_shuffle(&mut bucket);
                    let take = needed.min(bucket.len());
                    for compound in bucket.into_iter().take(take) {
                        picked_sequences.insert(compound.sequence.clone());
                        picked.push(compound);
                    }
                }
            }
        }
    }

    // 4) assemble the output experiment
    let mut out = LightTargetedExperiment::default();
    let pep_ids: HashSet<&str> = picked.iter().map(|c| c.id.as_str()).collect();
    for tr in &exp.transitions {
        if pep_ids.contains(tr.peptide_ref.as_str()) && !tr.flags.decoy() {
            out.transitions.push(tr.clone());
        }
    }
    let prot_ids: HashSet<&str> = picked
        .iter()
        .flat_map(|c| c.protein_refs.iter().map(String::as_str))
        .collect();
    for protein in &exp.proteins {
        if prot_ids.contains(protein.id.as_str()) {
            out.proteins.push(protein.clone());
        }
    }
    out.compounds = picked;
    Ok(out)
}

/// A random 32-bit seed, the port's `std::random_device{}()`: the low word of
/// the standard library's randomly keyed hasher.
fn random_seed() -> u32 {
    use std::hash::{BuildHasher, Hasher};
    let mut hasher = std::collections::hash_map::RandomState::new().build_hasher();
    hasher.write_u8(0);
    hasher.finish() as u32
}

/// The retention time of the best-scoring feature of each transition group,
/// keyed by the compound reference of the group's first transition.
///
/// Groups without features are skipped. With `use_qual_cutoff`, a group whose
/// best feature has an overall quality below `qual_cutoff` is skipped too. The
/// best feature is [`MRMTransitionGroup::best_feature`] (highest overall
/// quality, first on ties). When two groups share a compound reference, the
/// later group in key order wins.
///
/// # Errors
///
/// [`Error::MissingInformation`] when a group with features has no transition
/// (the source reads `getTransitions()[0]` unchecked).
pub fn simple_find_best_feature<C>(
    transition_group_map: &BTreeMap<String, MRMTransitionGroup<C, LightTransition>>,
    use_qual_cutoff: bool,
    qual_cutoff: f64,
) -> Result<BTreeMap<String, f64>> {
    let mut result = BTreeMap::new();
    for (key, group) in transition_group_map {
        if group.features().is_empty() {
            continue;
        }
        let best = group.best_feature()?;
        if use_qual_cutoff && f64::from(best.feature.base.quality) < qual_cutoff {
            continue;
        }
        let first = group.transitions().first().ok_or_else(|| {
            Error::MissingInformation(format!(
                "transition group '{key}' has features but no transitions"
            ))
        })?;
        result.insert(first.peptide_ref.clone(), best.rt);
    }
    Ok(result)
}
