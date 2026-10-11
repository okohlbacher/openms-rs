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
//! Not ported yet:
//!
//! - the `OpenMS::TargetedExperiment` overload of `selectSwathTransitions` and
//!   the corresponding instantiation of `checkSwathMapAndSelectTransitions`:
//!   the heavyweight `ANALYSIS/TARGETED/TargetedExperiment` is not in the crate;
//! - `sampleExperiment`: it draws through `Math::RandomShuffler`, whose
//!   source-exact engine is private to
//!   [`chemistry`](crate::chemistry) today, and seeds from
//!   `std::random_device` when the seed is 0.

use super::swath_map::SwathMap;
use super::transition_experiment::{LightTargetedExperiment, LightTransition};
use crate::kernel::MSExperiment;
use crate::kernel::mrm::MRMTransitionGroup;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

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
pub fn check_swath_map_and_select_transitions(
    exp: &MSExperiment,
    targeted_exp: &LightTargetedExperiment,
    selected: &mut LightTargetedExperiment,
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
    select_swath_transitions(targeted_exp, selected, min_upper_edge_dist, lower, upper);
    Ok(!selected.transitions.is_empty())
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
