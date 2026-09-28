// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Biosaur2 isotope pattern assembly, the source `detectIsotopePatterns_` and
//! its stages: the diagnostic initial calibration
//! (`performInitialIsotopeCalibration_`), the binned m/z lookup
//! (`buildFastMzLookup_`), candidate generation (`generateIsotopeCandidates_`),
//! RT-apex gating (`applyRtFiltering_`), mass-error calibration
//! (`refineIsotopeCalibration_`, `filterByCalibration_`) and the greedy
//! selection of non-overlapping patterns (`selectNonOverlappingPatterns_`),
//! with the averagine model (`computeAveragine_`) and the two cosine
//! correlations they rely on.
//!
//! # The averagine table
//!
//! The source builds a static table of C-only binomial isotope probabilities
//! with Boost.Math's `binomial_distribution` (`n_C = round(m / 111.1254 *
//! 4.9384)` carbons, at least one, for the neutral masses 0, 100, ..., 19900 Da,
//! `p = 0.0107`, isotopes 0-9, normalised). Boost evaluates the probability
//! mass through its incomplete-beta derivative with a Lanczos approximation and
//! glibc's `pow`, `exp` and `log`, so no crate reproduces its bits; `statrs`,
//! the crate the port already carries for Boost.Math, evaluates
//! `exp(ln C(n, k) + k ln p + (n - k) ln(1 - p))` instead. The table is
//! therefore the Release build's own, dumped by the oracle and embedded
//! (`averagine_table.rs`). `tests/biosaur2.rs` checks that it scales to
//! exactly what the library's `computeAveragine_` returns on every bin, and
//! that it agrees with an independent `statrs` evaluation to 1e-12.

use super::averagine_table::AVERAGINE_TABLE;
use super::{
    Biosaur2Algorithm, Biosaur2Message, Hill, IsotopeCandidate, PeptideFeature, g6, info, warn,
};
use crate::concept::constants::{C13C12_MASSDIFF_U, PROTON_MASS_U};
use crate::math::source_sort::source_sort_by;
use crate::processing::noise_estimation::x86::cvttsd2si32;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// `std::max(a, b)` on doubles: `b` only when `a < b`.
fn std_max(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// Per-isotope `(shift, sigma)` calibration in ppm, isotopes 1 to 9.
pub type IsotopeCalibration = BTreeMap<i32, (f64, f64)>;

/// Lightweight index entry for fast m/z-based hill lookup, the source private
/// `FastHillEntry`: the hill's position in the hill vector with its first and
/// last scan index, for quick RT overlap checks.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FastHillEntry {
    /// Position of the hill in the hill vector.
    pub hill_index: usize,
    /// First scan index of the hill.
    pub first_scan: usize,
    /// Last scan index of the hill.
    pub last_scan: usize,
}

/// Internal representation of a candidate isotope pattern, the source private
/// `PatternCandidate`: a monoisotopic hill with a set of isotope candidates,
/// the charge state and quality metrics used during pattern refinement.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PatternCandidate {
    /// Position of the monoisotopic hill in the hill vector.
    pub mono_index: usize,
    /// Monoisotopic m/z.
    pub mono_mz: f64,
    /// Charge state of the pattern.
    pub charge: i32,
    /// Cosine correlation in isotope-intensity space.
    pub cos_cor_isotopes: f64,
    /// Associated isotope candidates, in isotope order.
    pub isotopes: Vec<IsotopeCandidate>,
    /// Number of scans contributing to the monoisotopic hill.
    pub n_scans: usize,
}

/// Averagine-based theoretical isotope intensities, the source
/// `computeAveragine_`.
///
/// Uses the C-only binomial model on the 100 Da neutral-mass grid (the bin
/// `floor(neutral_mass / 100)`, clamped to 0-199; a mass whose bin does not fit
/// an `int` lands in bin 0, as `cvttsd2si` gives the Release build) and rescales
/// the probabilities so that the monoisotopic one equals `apex_intensity`.
/// Returns the ten intensities and the index of the expected most intense
/// isotope (at least 4).
///
/// The callers pass `mono_mz * charge`, which is not the neutral mass: the
/// proton masses are not subtracted, as in the source.
pub fn compute_averagine(neutral_mass: f64, apex_intensity: f64) -> ([f64; 10], usize) {
    let bin = cvttsd2si32((neutral_mass / 100.0).floor()).clamp(0, 199) as usize;
    let (_, max_pos, bits) = AVERAGINE_TABLE[bin];
    let probs = bits.map(f64::from_bits);
    let mono_prob = if probs[0] <= 0.0 { 1.0 } else { probs[0] };
    (probs.map(|p| apex_intensity * p / mono_prob), max_pos)
}

/// The number of carbons and the normalised probabilities of one averagine
/// bin, for the independent check of the embedded table.
pub fn averagine_bin(bin: usize) -> Option<(u32, [f64; 10])> {
    AVERAGINE_TABLE
        .get(bin)
        .map(|(n_c, _, bits)| (*n_c, bits.map(f64::from_bits)))
}

/// Cosine correlation between two 1-D intensity vectors over their common
/// length, the source `cosineCorrelation1D_`; 0 for an empty overlap or a
/// zero norm.
///
/// The source loop carries `#pragma omp simd reduction(+:dot,n1,n2)`, and the
/// Release build vectorises it over pairs of doubles: even and odd positions
/// accumulate separately, an odd last element joins the even lane, and each
/// sum is `(even + 0.0) + odd` (`libOpenMS.so` `0x19134d0`). This sums in that
/// order, so the value is the Release build's rather than the one a serial loop
/// gives.
pub fn cosine_correlation_1d(v1: &[f64], v2: &[f64]) -> f64 {
    let n = v1.len().min(v2.len());
    if n == 0 {
        return 0.0;
    }
    let mut dot = [0.0f64; 2];
    let mut n1 = [0.0f64; 2];
    let mut n2 = [0.0f64; 2];
    for pair in 0..n / 2 {
        for lane in 0..2 {
            let (a, b) = (v1[2 * pair + lane], v2[2 * pair + lane]);
            dot[lane] += a * b;
            n1[lane] += a * a;
            n2[lane] += b * b;
        }
    }
    if n % 2 == 1 {
        let (a, b) = (v1[n - 1], v2[n - 1]);
        dot[0] = a * b + dot[0];
        n1[0] = a * a + n1[0];
        n2[0] = b * b + n2[0];
    }
    let dot = (dot[0] + 0.0) + dot[1];
    let n1 = (n1[0] + 0.0) + n1[1];
    let n2 = (n2[0] + 0.0) + n2[1];
    if n1 == 0.0 || n2 == 0.0 {
        return 0.0;
    }
    dot / (n1.sqrt() * n2.sqrt())
}

/// Best cosine correlation of an experimental isotope envelope with the
/// averagine envelope and its truncation, the source
/// `checkingCosCorrelationForCarbon_` (the Python
/// `checking_cos_correlation_for_carbon`).
///
/// Truncations of `exp_full` from its full length down to two are tried in
/// turn; the first whose averagine prefix explains at least half of the whole
/// averagine intensity and whose correlation reaches `thresh` is returned with
/// its length. `(0, 1)` is returned for fewer than two experimental values, an
/// empty or non-positive averagine, or when no truncation qualifies.
pub fn checking_cos_correlation_for_carbon(
    theor_full: &[f64],
    exp_full: &[f64],
    thresh: f64,
) -> (f64, usize) {
    if exp_full.len() <= 1 || theor_full.is_empty() {
        return (0.0, 1);
    }
    let theor_total_sum = theor_full.iter().fold(0.0, |acc, &v| acc + v);
    if theor_total_sum <= 0.0 {
        return (0.0, 1);
    }
    let mut best_cor = 0.0;
    let mut best_pos = 1usize;
    let mut pos = exp_full.len();
    while pos != 1 {
        let suit_len = pos.min(theor_full.len());
        let theor = &theor_full[..suit_len];
        let exp = &exp_full[..pos];
        let theor_partial_sum = theor.iter().fold(0.0, |acc, &v| acc + v);
        let averagine_explained = if theor_total_sum > 0.0 {
            theor_partial_sum / theor_total_sum
        } else {
            0.0
        };
        let cor = cosine_correlation_1d(theor, exp);
        if averagine_explained >= 0.5 && cor >= thresh {
            if cor > best_cor {
                best_cor = cor;
                best_pos = pos;
            }
            break;
        }
        pos -= 1;
    }
    (best_cor, best_pos)
}

/// Whether a scan list is strictly ascending, so that the map of the source's
/// cosine is the list itself.
fn strictly_ascending(scans: &[usize]) -> bool {
    scans.windows(2).all(|w| w[0] < w[1])
}

/// Cosine correlation between two intensity traces over the scans they share,
/// the source `cosineCorrelation_`: the dot product runs over common scan
/// indices only, the norms over each whole trace; 0 for a zero norm.
///
/// The source keys each trace by scan in a `std::map` (a repeated scan keeps
/// its last intensity) and sums in ascending scan order. A hill's scans are
/// strictly ascending, and then the traces are summed directly in that order,
/// which gives the same value without building the maps.
pub fn cosine_correlation(
    intensities1: &[f64],
    scans1: &[usize],
    intensities2: &[f64],
    scans2: &[usize],
) -> f64 {
    let (mut dot, mut norm1, mut norm2) = (0.0, 0.0, 0.0);
    let n1 = scans1.len().min(intensities1.len());
    let n2 = scans2.len().min(intensities2.len());
    if strictly_ascending(&scans1[..n1]) && strictly_ascending(&scans2[..n2]) {
        let mut j = 0;
        for i in 0..n1 {
            while j < n2 && scans2[j] < scans1[i] {
                j += 1;
            }
            if j < n2 && scans2[j] == scans1[i] {
                dot += intensities1[i] * intensities2[j];
            }
            norm1 += intensities1[i] * intensities1[i];
        }
        for &value in &intensities2[..n2] {
            norm2 += value * value;
        }
    } else {
        let map1: BTreeMap<usize, f64> = scans1
            .iter()
            .copied()
            .zip(intensities1.iter().copied())
            .collect();
        let map2: BTreeMap<usize, f64> = scans2
            .iter()
            .copied()
            .zip(intensities2.iter().copied())
            .collect();
        for (scan, &i1) in &map1 {
            if let Some(&i2) = map2.get(scan) {
                dot += i1 * i2;
            }
            norm1 += i1 * i1;
        }
        for &i2 in map2.values() {
            norm2 += i2 * i2;
        }
    }
    if norm1 == 0.0 || norm2 == 0.0 {
        return 0.0;
    }
    dot / (norm1.sqrt() * norm2.sqrt())
}

/// Recommended number of isotopes to keep before a valley in the isotope
/// envelope, the source `checkIsotopeValleySplit_`.
///
/// The apex intensities of the isotope hills are smoothed ([`mean_filter`],
/// half-width 3); from position `max(4, argmax + 1)` on, the first position
/// whose value times `ivf` is below the maximum to its right is returned, or
/// the number of isotopes when there is none. One isotope or none is returned
/// unchanged.
///
/// The source defines this helper but never calls it: its candidate
/// generation applies the isotope valley factor inline. It is ported for
/// completeness.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when an isotope names a hill that `hills`
/// does not hold; the source dereferences the end iterator.
///
/// [`mean_filter`]: super::hills::mean_filter
pub fn check_isotope_valley_split(
    isotopes: &[IsotopeCandidate],
    hills: &[Hill],
    ivf: f64,
) -> Result<usize> {
    if isotopes.len() <= 1 {
        return Ok(isotopes.len());
    }
    let mut intensities = Vec::with_capacity(isotopes.len());
    for iso in isotopes {
        let hill = hills
            .iter()
            .find(|h| h.hill_idx == iso.hill_idx)
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "isotope hill {} is not among the hills",
                    iso.hill_idx
                ))
            })?;
        intensities.push(hill.intensity_apex);
    }
    let smoothed = super::hills::mean_filter(&intensities, 3);
    let max_pos = first_max(&smoothed);
    let min_check_pos = 4usize.max(max_pos + 1);
    for i in min_check_pos..smoothed.len() - 1 {
        let local_min = smoothed[i];
        let right_max = smoothed[first_max(&smoothed[i + 1..]) + i + 1];
        if local_min * ivf < right_max {
            return Ok(i);
        }
    }
    Ok(smoothed.len())
}

/// `std::max_element`: the position of the first maximum under `<`.
fn first_max(values: &[f64]) -> usize {
    let mut best = 0;
    for (i, &value) in values.iter().enumerate() {
        if values[best] < value {
            best = i;
        }
    }
    best
}

/// The initial isotope calibration on the raw hills, the source
/// `performInitialIsotopeCalibration_`.
///
/// Scans all hills (sorted by m/z) for regularly spaced C13 isotope peaks,
/// from the highest charge down until a first isotope is found, and derives
/// per-isotope shift/sigma estimates in ppm: [`calibrate_mass`] for isotopes 1
/// to 9 with at least 1000 errors, an extrapolation from the two previous
/// isotopes for the others from 4 on, `(0, itol)` otherwise. The source
/// computes this for its log lines only and discards the result, as the Python
/// reference does.
///
/// # Errors
///
/// As [`calibrate_mass`].
///
/// [`calibrate_mass`]: super::hills::calibrate_mass
pub fn perform_initial_isotope_calibration(
    hills: &[Hill],
    itol_ppm: f64,
    min_charge: i32,
    max_charge: i32,
    enable_isotope_calib: bool,
    messages: &mut Vec<Biosaur2Message>,
) -> Result<IsotopeCalibration> {
    let mut calib: IsotopeCalibration = (1..=9).map(|ic| (ic, (0.0, itol_ppm))).collect();
    if !enable_isotope_calib || hills.is_empty() {
        return Ok(calib);
    }
    info(messages, "Performing isotope calibration...".to_owned());
    let mut errors: BTreeMap<i32, Vec<f64>> = (1..=9).map(|ic| (ic, Vec::new())).collect();
    for (i, mono_hill) in hills.iter().enumerate() {
        let mono_mz = mono_hill.mz_weighted_mean;
        let mut charge = max_charge;
        while charge >= min_charge {
            let mz_spacing = C13C12_MASSDIFF_U / f64::from(charge);
            let mut found_first = false;
            for iso_num in 1..=9i32 {
                let expected_mz = mono_mz + f64::from(iso_num) * mz_spacing;
                let mz_tolerance = expected_mz * itol_ppm * 1e-6;
                let mut best: Option<usize> = None;
                let mut best_intensity = -1.0;
                for (j, hill) in hills.iter().enumerate().skip(i + 1) {
                    if hill.mz_weighted_mean > expected_mz + mz_tolerance {
                        break;
                    }
                    if hill.mz_weighted_mean < expected_mz - mz_tolerance {
                        continue;
                    }
                    let diff = (hill.mz_weighted_mean - expected_mz).abs();
                    if diff <= mz_tolerance && hill.intensity_apex > best_intensity {
                        best_intensity = hill.intensity_apex;
                        best = Some(j);
                    }
                }
                if let Some(j) = best {
                    if mono_hill.length >= 3 {
                        let ppm =
                            super::hills::calculate_ppm(hills[j].mz_weighted_mean, expected_mz);
                        errors.entry(iso_num).or_default().push(ppm);
                        if iso_num == 1 {
                            found_first = true;
                        }
                    }
                }
            }
            if found_first {
                break;
            }
            charge -= 1;
        }
    }
    for ic in 1..=3 {
        if errors[&ic].len() >= 1000 {
            let (shift, sigma) = super::hills::calibrate_mass(&errors[&ic], 0.05)?;
            calib.insert(ic, (shift, sigma));
            info(
                messages,
                format!(
                    "Isotope {ic} calibration: shift={} ppm, sigma={} ppm",
                    g6(shift),
                    g6(sigma)
                ),
            );
        }
    }
    extrapolate(&mut calib, &errors, itol_ppm)?;
    let (shift, sigma) = calib[&1];
    info(
        messages,
        format!(
            "Isotope 1 calibration: shift={} ppm, sigma={} ppm",
            g6(shift),
            g6(sigma)
        ),
    );
    Ok(calib)
}

/// Isotopes 4 to 9: their own calibration with at least 1000 errors, else the
/// previous isotope's shift plus its step from the one before, and its sigma
/// times its ratio to the one before (floored at 0.1).
fn extrapolate(
    calib: &mut IsotopeCalibration,
    errors: &BTreeMap<i32, Vec<f64>>,
    itol_ppm: f64,
) -> Result<()> {
    for ic in 4..=9 {
        if errors[&ic].len() >= 1000 {
            let fit = super::hills::calibrate_mass(&errors[&ic], 0.05)?;
            calib.insert(ic, fit);
        } else {
            let prev = calib[&(ic - 1)];
            let prev2 = calib.get(&(ic - 2)).copied().unwrap_or((0.0, itol_ppm));
            let shift_delta = prev.0 - prev2.0;
            let sigma_ratio = prev.1 / std_max(prev2.1, 0.1);
            calib.insert(ic, (prev.0 + shift_delta, prev.1 * sigma_ratio));
        }
    }
    Ok(())
}

/// The fast m/z lookup and the ion-mobility bins of the hills, the source
/// `buildFastMzLookup_`.
///
/// Every hill with scans is registered in its m/z bin and both neighbours
/// (bins `htol * 1e-6 * max(mz)` wide, the uncalibrated `htol`); with `use_im`
/// a hill with a non-negative ion mobility gets the bin
/// `ion_mobility / paseftol`, and -1 (no bin) otherwise. Returns the bin width
/// (0 when no hill has a positive m/z, which puts every hill in bin 0), the
/// lookup and the bins.
pub fn build_fast_mz_lookup(
    algorithm: &Biosaur2Algorithm,
    hills: &[Hill],
    use_im: bool,
) -> (f64, BTreeMap<i32, Vec<FastHillEntry>>, Vec<i32>) {
    let mut lookup: BTreeMap<i32, Vec<FastHillEntry>> = BTreeMap::new();
    let mut im_bins = vec![-1i32; hills.len()];
    let mut max_mz_value = 0.0;
    for hill in hills {
        max_mz_value = std_max(max_mz_value, hill.mz_weighted_mean);
    }
    let mz_step = if max_mz_value > 0.0 {
        algorithm.htol * 1e-6 * max_mz_value
    } else {
        0.0
    };
    for (index, hill) in hills.iter().enumerate() {
        let (Some(&first_scan), Some(&last_scan)) =
            (hill.scan_indices.first(), hill.scan_indices.last())
        else {
            continue;
        };
        let bin = if mz_step > 0.0 {
            cvttsd2si32(hill.mz_weighted_mean / mz_step)
        } else {
            0
        };
        let entry = FastHillEntry {
            hill_index: index,
            first_scan,
            last_scan,
        };
        for key in [bin.wrapping_sub(1), bin, bin.wrapping_add(1)] {
            lookup.entry(key).or_default().push(entry);
        }
        if use_im && hill.ion_mobility_median >= 0.0 {
            im_bins[index] = cvttsd2si32(hill.ion_mobility_median / algorithm.paseftol);
        }
    }
    (mz_step, lookup, im_bins)
}

/// Whether two ascending scan lists share a scan.
fn has_scan_overlap(a: &[usize], b: &[usize]) -> bool {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i] == b[j] {
            return true;
        }
        if a[i] < b[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    false
}

/// The charges a pattern of each charge bans for lower charges of the same
/// monoisotopic hill, the source's `charge_ban_map` (from the Python
/// reference).
fn charge_bans(charge: i32) -> &'static [i32] {
    match charge {
        8 => &[1, 2, 4],
        7 | 5 | 3 | 2 | 1 => &[1],
        6 => &[1, 2, 3],
        4 => &[1, 2],
        _ => &[],
    }
}

/// Generate the initial isotope pattern candidates for every monoisotopic hill,
/// the source `generateIsotopeCandidates_` (the Python `get_initial_isotopes`).
///
/// For each hill and each charge from `max_charge` down to `min_charge`,
/// isotopes 1 to 9 are looked up in the fast m/z lookup: a candidate is another
/// hill in the expected isotope's bin, overlapping in scans (and, with
/// `use_im`, at most one ion-mobility bin away), within `itol_ppm` of the
/// expected m/z and with an RT-profile cosine of at least 0.6. The search stops
/// at the first isotope without a candidate. Every combination of one candidate
/// per isotope is scored against the averagine envelope: from the expected most
/// intense isotope on, the envelope is cut before an isotope that rises to
/// `ivf` times the lowest one since; a combination whose cosine reaches 0.6 on a
/// prefix explaining half of the averagine becomes a pattern of that length.
/// A pattern of one charge requires patterns of the charges it bans (see the
/// source's charge ban map) to reach its length.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when one monoisotopic hill and charge would
/// enumerate more than [`Biosaur2Limits::max_isotope_combinations`]
/// combinations, checked before the enumeration.
///
/// [`Biosaur2Limits::max_isotope_combinations`]: super::Biosaur2Limits::max_isotope_combinations
#[allow(clippy::too_many_arguments)]
pub fn generate_isotope_candidates(
    algorithm: &Biosaur2Algorithm,
    hills: &[Hill],
    itol_ppm: f64,
    min_charge: i32,
    max_charge: i32,
    ivf: f64,
    mz_step: f64,
    hills_mz_fast: &BTreeMap<i32, Vec<FastHillEntry>>,
    hill_idx_to_index: &BTreeMap<usize, usize>,
    hill_im_bins: &[i32],
    use_im: bool,
) -> Result<Vec<PatternCandidate>> {
    let mut ready = Vec::new();
    if hills.is_empty() {
        return Ok(ready);
    }
    let charges: Vec<i32> = if min_charge <= max_charge {
        (min_charge..=max_charge).rev().collect()
    } else {
        Vec::new()
    };
    for (i, mono_hill) in hills.iter().enumerate() {
        let (Some(&first1), Some(&last1)) = (
            mono_hill.scan_indices.first(),
            mono_hill.scan_indices.last(),
        ) else {
            continue;
        };
        let mono_mz = mono_hill.mz_weighted_mean;
        let mz_tol = itol_ppm * 1e-6 * mono_mz;
        let mut banned_charges: BTreeMap<i32, usize> = BTreeMap::new();
        for &charge in &charges {
            let mut candidates: Vec<Vec<IsotopeCandidate>> = Vec::new();
            for iso_num in 1..=9usize {
                let mut tmp = Vec::new();
                let expected_mz = mono_mz + C13C12_MASSDIFF_U * iso_num as f64 / f64::from(charge);
                let fast = if mz_step > 0.0 {
                    cvttsd2si32(expected_mz / mz_step)
                } else {
                    0
                };
                if let Some(entries) = hills_mz_fast.get(&fast) {
                    for entry in entries {
                        let idx2 = entry.hill_index;
                        if idx2 == i {
                            continue;
                        }
                        let hill2 = &hills[idx2];
                        if hill2.scan_indices.is_empty() {
                            continue;
                        }
                        if last1 < entry.first_scan || entry.last_scan < first1 {
                            continue;
                        }
                        if use_im
                            && hill_im_bins.len() > i
                            && hill_im_bins.len() > idx2
                            && hill_im_bins[i] >= 0
                            && hill_im_bins[idx2] >= 0
                            && (i64::from(hill_im_bins[i]) - i64::from(hill_im_bins[idx2])).abs()
                                > 1
                        {
                            continue;
                        }
                        let mass_diff_abs = hill2.mz_weighted_mean - expected_mz;
                        if mass_diff_abs.abs() > mz_tol {
                            continue;
                        }
                        if !has_scan_overlap(&mono_hill.scan_indices, &hill2.scan_indices) {
                            continue;
                        }
                        let cos_cor_rt = cosine_correlation(
                            &mono_hill.intensities,
                            &mono_hill.scan_indices,
                            &hill2.intensities,
                            &hill2.scan_indices,
                        );
                        if cos_cor_rt < 0.6 {
                            continue;
                        }
                        tmp.push(IsotopeCandidate {
                            hill_idx: hill2.hill_idx,
                            isotope_number: iso_num,
                            mass_diff_ppm: mass_diff_abs * 1e6 / expected_mz,
                            cos_corr: cos_cor_rt,
                        });
                    }
                }
                if !tmp.is_empty() {
                    candidates.push(tmp);
                }
                if candidates.len() < iso_num {
                    break;
                }
            }

            let min_required = banned_charges.get(&charge).copied().unwrap_or(1);
            if candidates.len() < min_required || candidates.is_empty() {
                continue;
            }
            let mut combinations = 1usize;
            for options in &candidates {
                combinations = combinations.saturating_mul(options.len());
            }
            if combinations > algorithm.limits.max_isotope_combinations {
                return Err(Error::InvalidValue(format!(
                    "Biosaur2: {combinations} isotope combinations for one hill and charge exceed the limit of {}",
                    algorithm.limits.max_isotope_combinations
                )));
            }
            let apex1 = mono_hill.intensity_apex;
            let (theoretical, max_pos) = compute_averagine(mono_mz * f64::from(charge), apex1);

            let mut indices = vec![0usize; candidates.len()];
            loop {
                let picked: Vec<IsotopeCandidate> = candidates
                    .iter()
                    .zip(&indices)
                    .map(|(options, &k)| options[k])
                    .collect();
                let mut exp_intensity = Vec::with_capacity(picked.len() + 1);
                exp_intensity.push(apex1);
                let mut local_minimum = 0.0;
                let mut local_minimum_pos = 0usize;
                let mut i_local_isotope = 1usize;
                for iso in &picked {
                    let apex2 = hill_idx_to_index
                        .get(&iso.hill_idx)
                        .map_or(0.0, |&position| hills[position].intensity_apex);
                    if i_local_isotope > max_pos {
                        if i_local_isotope == max_pos + 1 || apex2 < local_minimum {
                            local_minimum = apex2;
                            local_minimum_pos = i_local_isotope;
                        }
                        if apex2 >= ivf * local_minimum {
                            if local_minimum_pos + 1 < exp_intensity.len() {
                                exp_intensity.truncate(local_minimum_pos + 1);
                            }
                            break;
                        }
                    }
                    exp_intensity.push(apex2);
                    i_local_isotope += 1;
                }
                let (cos_corr_iso, best_pos) =
                    checking_cos_correlation_for_carbon(&theoretical, &exp_intensity, 0.6);
                if cos_corr_iso > 0.0 && best_pos > 1 {
                    let n_iso_used = (best_pos - 1).min(picked.len());
                    let pattern = PatternCandidate {
                        mono_index: i,
                        mono_mz,
                        charge,
                        cos_cor_isotopes: cos_corr_iso,
                        isotopes: picked[..n_iso_used].to_vec(),
                        n_scans: mono_hill.length,
                    };
                    if !pattern.isotopes.is_empty() {
                        ready.push(pattern);
                        for &banned in charge_bans(charge) {
                            let slot = banned_charges.entry(banned).or_insert(0);
                            *slot = (*slot).max(n_iso_used);
                        }
                    }
                }
                // Next combination: the first isotope order varies fastest.
                let mut position = 0;
                while position < indices.len() {
                    indices[position] += 1;
                    if indices[position] < candidates[position].len() {
                        break;
                    }
                    indices[position] = 0;
                    position += 1;
                }
                if position == indices.len() {
                    break;
                }
            }
        }
    }
    Ok(ready)
}

/// RT-apex gating of the pattern candidates, the source `applyRtFiltering_`:
/// isotopes whose hill apex lies more than `hrttol` seconds from the
/// monoisotopic apex are dropped, and patterns left without isotopes with them.
/// The candidates are returned unchanged when `hrttol` is 0 or there are none.
pub fn apply_rt_filtering(
    algorithm: &Biosaur2Algorithm,
    candidates: &[PatternCandidate],
    hills: &[Hill],
    hill_idx_to_index: &BTreeMap<usize, usize>,
    messages: &mut Vec<Biosaur2Message>,
) -> Vec<PatternCandidate> {
    if algorithm.hrttol <= 0.0 || candidates.is_empty() {
        return candidates.to_vec();
    }
    let mut filtered = Vec::with_capacity(candidates.len());
    for candidate in candidates {
        let mono_rt_apex = hills[candidate.mono_index].rt_apex;
        let kept: Vec<IsotopeCandidate> = candidate
            .isotopes
            .iter()
            .filter(|iso| {
                hill_idx_to_index
                    .get(&iso.hill_idx)
                    .is_some_and(|&position| {
                        (hills[position].rt_apex - mono_rt_apex).abs() <= algorithm.hrttol
                    })
            })
            .copied()
            .collect();
        if !kept.is_empty() {
            filtered.push(PatternCandidate {
                isotopes: kept,
                ..candidate.clone()
            });
        }
    }
    info(
        messages,
        format!(
            "After RT apex filter (hrttol={} s), {} potential isotope patterns remain.",
            g6(algorithm.hrttol),
            filtered.len()
        ),
    );
    filtered
}

/// The isotope mass-error calibration from the pattern candidates, the source
/// `refineIsotopeCalibration_`: the mass errors of the isotopes of every
/// candidate whose monoisotopic hill spans at least 3 scans, fitted with
/// [`calibrate_mass`](super::hills::calibrate_mass) per isotope with at least
/// 1000 errors and extrapolated for isotopes 4 to 9 otherwise; `(0, itol)`
/// everywhere when the calibration is disabled.
///
/// # Errors
///
/// As [`calibrate_mass`](super::hills::calibrate_mass).
pub fn refine_isotope_calibration(
    candidates: &[PatternCandidate],
    itol_ppm: f64,
    enable_isotope_calib: bool,
) -> Result<IsotopeCalibration> {
    let mut errors: BTreeMap<i32, Vec<f64>> = (1..=9).map(|ic| (ic, Vec::new())).collect();
    for candidate in candidates {
        if candidate.n_scans < 3 {
            continue;
        }
        for iso in &candidate.isotopes {
            if (1..=9).contains(&iso.isotope_number) {
                errors
                    .entry(iso.isotope_number as i32)
                    .or_default()
                    .push(iso.mass_diff_ppm);
            }
        }
    }
    let mut calib: IsotopeCalibration = (1..=9).map(|ic| (ic, (0.0, itol_ppm))).collect();
    if !enable_isotope_calib {
        return Ok(calib);
    }
    for ic in 1..=3 {
        if errors[&ic].len() >= 1000 {
            calib.insert(ic, super::hills::calibrate_mass(&errors[&ic], 0.05)?);
        }
    }
    extrapolate(&mut calib, &errors, itol_ppm)?;
    Ok(calib)
}

/// Filter the pattern candidates with the isotope calibration and re-check
/// their averagine correlation, the source `filterByCalibration_`.
///
/// With calibration enabled each pattern keeps its leading isotopes whose mass
/// error lies within 5 sigma of the isotope's shift. The averagine cosine is
/// then recomputed on the remaining envelope ([`checking_cos_correlation_for_carbon`])
/// and the pattern truncated to its best length; patterns without isotopes or
/// without a qualifying truncation are dropped.
pub fn filter_by_calibration(
    candidates: &[PatternCandidate],
    hills: &[Hill],
    hill_idx_to_index: &BTreeMap<usize, usize>,
    calibration: &IsotopeCalibration,
    enable_isotope_calib: bool,
) -> Vec<PatternCandidate> {
    let mut filtered = Vec::with_capacity(candidates.len());
    for candidate in candidates {
        let mut pattern = candidate.clone();
        if enable_isotope_calib {
            let mut kept = Vec::new();
            for iso in &pattern.isotopes {
                let Some(&(shift, sigma)) = i32::try_from(iso.isotope_number)
                    .ok()
                    .and_then(|number| calibration.get(&number))
                else {
                    continue;
                };
                if (iso.mass_diff_ppm - shift).abs() <= 5.0 * sigma {
                    kept.push(*iso);
                } else {
                    break;
                }
            }
            pattern.isotopes = kept;
        }
        if pattern.isotopes.is_empty() {
            continue;
        }
        let mono_hill = &hills[pattern.mono_index];
        let apex1 = mono_hill.intensity_apex;
        let (theoretical, _) = compute_averagine(
            mono_hill.mz_weighted_mean * f64::from(pattern.charge),
            apex1,
        );
        let mut exp_intensity = Vec::with_capacity(pattern.isotopes.len() + 1);
        exp_intensity.push(apex1);
        for iso in &pattern.isotopes {
            exp_intensity.push(
                hill_idx_to_index
                    .get(&iso.hill_idx)
                    .map_or(0.0, |&position| hills[position].intensity_apex),
            );
        }
        let (cos_corr_iso, best_pos) =
            checking_cos_correlation_for_carbon(&theoretical, &exp_intensity, 0.6);
        if cos_corr_iso <= 0.0 || best_pos <= 1 {
            continue;
        }
        let n_iso_used = (best_pos - 1).min(pattern.isotopes.len());
        pattern.isotopes.truncate(n_iso_used);
        pattern.cos_cor_isotopes = cos_corr_iso;
        if !pattern.isotopes.is_empty() {
            filtered.push(pattern);
        }
    }
    filtered
}

/// Log obviously inconsistent isotope assignments, the source
/// `debugCheckIsotopeConsistency_`: an isotope hill whose RT apex lies more
/// than `hrttol` (+1e-6) seconds from the monoisotopic apex, or whose m/z lies
/// more than `max(80, 10 * itol)` ppm from the expected isotope position. Both
/// are warnings (`OPENMS_LOG_WARN`) labelled with `stage_label`.
#[allow(clippy::too_many_arguments)]
pub fn debug_check_isotope_consistency(
    algorithm: &Biosaur2Algorithm,
    stage_label: &str,
    mono_mz_center: f64,
    mono_rt_apex: f64,
    mono_hill_idx: usize,
    charge: i32,
    itol_ppm: f64,
    iso_hill: &Hill,
    isotope_number: usize,
    messages: &mut Vec<Biosaur2Message>,
) {
    if algorithm.hrttol > 0.0 {
        let rt_delta = (iso_hill.rt_apex - mono_rt_apex).abs();
        if rt_delta > algorithm.hrttol + 1e-6 {
            warn(
                messages,
                format!(
                    "Biosaur2 isotope debug ({stage_label}): mono m/z {} (charge {charge}, mono hill_idx={mono_hill_idx}) uses isotope hill_idx={} with RT apex delta {} s > hrttol={}. This indicates an inconsistent isotope assignment.",
                    g6(mono_mz_center),
                    iso_hill.hill_idx,
                    g6(rt_delta),
                    g6(algorithm.hrttol)
                ),
            );
        }
    }
    if charge > 0 {
        let expected_mz =
            mono_mz_center + isotope_number as f64 * C13C12_MASSDIFF_U / f64::from(charge);
        let observed_mz = iso_hill.mz_weighted_mean;
        let diff_ppm = (observed_mz - expected_mz) / expected_mz * 1e6;
        let threshold = std_max(80.0, 10.0 * itol_ppm);
        if diff_ppm.abs() > threshold {
            warn(
                messages,
                format!(
                    "Biosaur2 isotope debug ({stage_label}): mono m/z {} (charge {charge}, mono hill_idx={mono_hill_idx}) uses isotope #{isotope_number} (hill_idx={}) at m/z={} which is {} ppm away from expected {} (itol={} ppm).",
                    g6(mono_mz_center),
                    iso_hill.hill_idx,
                    g6(observed_mz),
                    g6(diff_ppm),
                    g6(expected_mz),
                    g6(itol_ppm)
                ),
            );
        }
    }
}

/// Greedily select non-overlapping isotope patterns and assemble the peptide
/// features, the source `selectNonOverlappingPatterns_` (the Python
/// `ready_final`).
///
/// The patterns are ranked by isotope count, then cosine, both descending (the
/// source's `std::sort`, ties as the Release build orders them). A pattern
/// whose monoisotopic hill is taken is skipped; one with a taken isotope hill
/// keeps its leading free isotopes, re-checked against the averagine envelope
/// and truncated. Each accepted pattern becomes a [`PeptideFeature`] (with the
/// apex and sum intensities of the first `iuse` isotopes added, all of them for
/// -1) and claims its hills.
///
/// # Errors
///
/// Only an introsort refusal of [`source_sort_by`], which these comparisons
/// cannot provoke.
pub fn select_non_overlapping_patterns(
    algorithm: &Biosaur2Algorithm,
    filtered_ready: &[PatternCandidate],
    hills: &[Hill],
    itol_ppm: f64,
    messages: &mut Vec<Biosaur2Message>,
) -> Result<Vec<PeptideFeature>> {
    let mut features = Vec::new();
    if filtered_ready.is_empty() {
        return Ok(features);
    }
    let mut sorted = filtered_ready.to_vec();
    source_sort_by(&mut sorted, |a, b| {
        if a.isotopes.len() != b.isotopes.len() {
            a.isotopes.len() > b.isotopes.len()
        } else {
            a.cos_cor_isotopes > b.cos_cor_isotopes
        }
    })?;
    let lookup: BTreeMap<usize, usize> = hills
        .iter()
        .enumerate()
        .map(|(position, hill)| (hill.hill_idx, position))
        .collect();
    let mut occupied: BTreeSet<usize> = BTreeSet::new();
    for candidate in &sorted {
        let mut pattern = candidate.clone();
        let mono_hill = &hills[pattern.mono_index];
        let mono_mz_center = mono_hill.mz_weighted_mean;
        if occupied.contains(&mono_hill.hill_idx) {
            continue;
        }
        if pattern
            .isotopes
            .iter()
            .any(|iso| occupied.contains(&iso.hill_idx))
        {
            let mut tmp_iso: Vec<IsotopeCandidate> = pattern
                .isotopes
                .iter()
                .take_while(|iso| !occupied.contains(&iso.hill_idx))
                .copied()
                .collect();
            if tmp_iso.is_empty() {
                continue;
            }
            let apex1 = mono_hill.intensity_apex;
            let (theor_full, _) = compute_averagine(
                mono_hill.mz_weighted_mean * f64::from(pattern.charge),
                apex1,
            );
            let len = (tmp_iso.len() + 1).min(theor_full.len());
            if len <= 1 {
                continue;
            }
            let theor = &theor_full[..len];
            let mut exp = vec![0.0; len];
            exp[0] = apex1;
            for k in 1..len {
                exp[k] = lookup
                    .get(&tmp_iso[k - 1].hill_idx)
                    .map_or(0.0, |&position| hills[position].intensity_apex);
            }
            let (cos_corr_iso, best_pos) = checking_cos_correlation_for_carbon(theor, &exp, 0.6);
            if cos_corr_iso <= 0.0 || best_pos <= 1 {
                continue;
            }
            let n_iso_used = (best_pos - 1).min(tmp_iso.len());
            tmp_iso.truncate(n_iso_used);
            pattern.isotopes = tmp_iso;
            pattern.cos_cor_isotopes = cos_corr_iso;
        }
        if pattern.isotopes.is_empty() {
            continue;
        }
        for iso in &pattern.isotopes {
            if let Some(&position) = lookup.get(&iso.hill_idx) {
                debug_check_isotope_consistency(
                    algorithm,
                    "detectIsotopePatterns_",
                    mono_mz_center,
                    mono_hill.rt_apex,
                    mono_hill.hill_idx,
                    pattern.charge,
                    itol_ppm,
                    &hills[position],
                    iso.isotope_number,
                    messages,
                );
            }
        }

        let mut feature = PeptideFeature {
            mz: pattern.mono_mz,
            rt_start: mono_hill.rt_start,
            rt_end: mono_hill.rt_end,
            rt_apex: mono_hill.rt_apex,
            intensity_apex: mono_hill.intensity_apex,
            intensity_sum: mono_hill.intensity_sum,
            ..PeptideFeature::default()
        };
        if algorithm.iuse != 0 {
            let available = pattern.isotopes.len();
            let to_add = if algorithm.iuse == -1 {
                available
            } else {
                available.min(usize::try_from(algorithm.iuse).unwrap_or(0))
            };
            for iso in &pattern.isotopes[..to_add] {
                if let Some(&position) = lookup.get(&iso.hill_idx) {
                    feature.intensity_apex += hills[position].intensity_apex;
                    feature.intensity_sum += hills[position].intensity_sum;
                }
            }
        }
        feature.charge = pattern.charge;
        feature.n_isotopes = pattern.isotopes.len() + 1;
        feature.n_scans = mono_hill.length;
        feature.isotopes = pattern.isotopes.clone();
        feature.mono_hill_idx = mono_hill.hill_idx;
        feature.drift_time = mono_hill.drift_time_median;
        feature.ion_mobility = mono_hill.ion_mobility_median;
        let charge = f64::from(pattern.charge);
        feature.mass_calib = if algorithm.negative_mode {
            pattern.mono_mz * charge + PROTON_MASS_U * charge
        } else {
            pattern.mono_mz * charge - PROTON_MASS_U * charge
        };
        features.push(feature);
        occupied.insert(mono_hill.hill_idx);
        for iso in &pattern.isotopes {
            occupied.insert(iso.hill_idx);
        }
    }
    Ok(features)
}

/// Detect isotope patterns and assemble peptide features, the source
/// `detectIsotopePatterns_`.
///
/// Sorts `hills` by weighted m/z (the source's `std::sort`, which is why the
/// hills of [`Biosaur2Output`](super::Biosaur2Output) come out in m/z order),
/// then runs the stages of this module in order: the diagnostic initial
/// calibration, the fast lookup, candidate generation, RT gating, the
/// calibration refinement and filter, and the greedy selection.
///
/// # Errors
///
/// As the stages.
pub fn detect_isotope_patterns(
    algorithm: &Biosaur2Algorithm,
    hills: &mut Vec<Hill>,
    use_im: bool,
    messages: &mut Vec<Biosaur2Message>,
) -> Result<Vec<PeptideFeature>> {
    info(messages, "Detecting isotope patterns...".to_owned());
    source_sort_by(hills, |a, b| a.mz_weighted_mean < b.mz_weighted_mean)?;
    let hill_idx_to_index: BTreeMap<usize, usize> = hills
        .iter()
        .enumerate()
        .map(|(position, hill)| (hill.hill_idx, position))
        .collect();
    let enable_isotope_calib = !algorithm.ignore_iso_calib;
    perform_initial_isotope_calibration(
        hills,
        algorithm.itol,
        algorithm.cmin,
        algorithm.cmax,
        enable_isotope_calib,
        messages,
    )?;
    let (mz_step, lookup, im_bins) = build_fast_mz_lookup(algorithm, hills, use_im);
    let ready = generate_isotope_candidates(
        algorithm,
        hills,
        algorithm.itol,
        algorithm.cmin,
        algorithm.cmax,
        algorithm.ivf,
        mz_step,
        &lookup,
        &hill_idx_to_index,
        &im_bins,
        use_im,
    )?;
    let rt_ready = apply_rt_filtering(algorithm, &ready, hills, &hill_idx_to_index, messages);
    drop(ready);
    let calibration = refine_isotope_calibration(&rt_ready, algorithm.itol, enable_isotope_calib)?;
    let filtered = filter_by_calibration(
        &rt_ready,
        hills,
        &hill_idx_to_index,
        &calibration,
        enable_isotope_calib,
    );
    let features =
        select_non_overlapping_patterns(algorithm, &filtered, hills, algorithm.itol, messages)?;
    info(
        messages,
        format!("Detected {} features with isotope patterns", features.len()),
    );
    Ok(features)
}
