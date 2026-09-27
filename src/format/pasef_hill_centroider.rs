// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ion-mobility-axis hill centroiding of one TIMS-PASEF frame, the source
//! namespace `Internal::PASEFHillCentroider`
//! (`FORMAT/HANDLERS/PASEFHillCentroider.h`).
//!
//! The algorithm is modelled on Biosaur2's `linkScanToHills_` and
//! `splitHills_`, with the IM scans of one frame playing the role that
//! consecutive RT spectra play in Biosaur2. [`centroid_frame`]
//!
//! 1. groups the peaks by IM scan (by scan index when given, otherwise by equal
//!    IM value),
//! 2. links peaks of consecutive IM scans by m/z proximity, assigning
//!    current-scan peaks greedily in descending intensity order,
//! 3. splits hills at intensity valleys,
//! 4. drops hills shorter than [`Params::min_hill_length`], and
//! 5. collapses each surviving hill to one [`Centroid`].
//!
//! The source is used by the Bruker TDF reader (`BrukerTimsFile`) for its
//! `HillBased` centroiding modes; it has no dependency on the reader and is
//! ported here on its own. See `docs/BRUKER_TIMS_SUPPORT.md` for the member
//! mapping and the evidence.
//!
//! # Native differences
//!
//! - Input is three (optionally four) slices instead of raw pointers plus a
//!   count; slices of different lengths are refused.
//! - Non-finite m/z, intensity or IM values, and non-finite parameters, are
//!   refused with [`Error::InvalidValue`]. The source sorts with `<`/`>`
//!   comparators, which is undefined behaviour for NaN, and casts
//!   `mz / mz_step` to `int`, which is undefined for a value out of range.
//! - Every sort is stable. The source uses `std::sort`, whose order among
//!   equal keys is unspecified; the port resolves such ties in input order.
//!   Ties only arise for peaks with identical sort keys (identical intensity
//!   within one scan, or identical `(scan, m/z)`), and the source makes no
//!   promise about them.
//! - A run of empty IM scans between two non-empty ones is aged in one step
//!   instead of scan by scan. This gives the same tails as the source (ageing
//!   and pruning are monotone and pruning keeps the tail order) while keeping
//!   the cost linear in the number of peaks; the source walks every scan index
//!   from the smallest to the largest, and with a largest index of
//!   `u32::MAX` its `s <= scan_max; ++s` loop never ends.

use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// Parameters of IM-axis hill-based centroiding of one TIMS-PASEF frame.
///
/// The source struct `PASEFHillCentroider::Params`; the defaults are the
/// source's.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Params {
    /// m/z linking tolerance in ppm. A peak in the next IM scan extends a hill
    /// when its m/z is within this tolerance of the previous-scan peak.
    /// Default 8.
    ///
    /// The m/z bin width of the linker is `mz_tol_ppm * 1e-6 * max_mz`, as in
    /// Biosaur2's `computeHillMzStep_`; a tolerance that makes it zero or
    /// negative yields no centroids, as in the source.
    pub mz_tol_ppm: f64,
    /// Hill-valley factor (hvf): a valley splits a hill only when both
    /// `left_max / valley` and `right_max / valley` reach this value. Smaller
    /// values split more aggressively. Default 1.3, Biosaur2's default.
    pub valley_factor: f64,
    /// Minimum number of peaks per hill; shorter hills are dropped. Default 1.
    ///
    /// The source's rationale for 1: detector centroiding produces sparse
    /// per-scan peaks and many real ions span only one IM scan of a frame
    /// (roughly 75% of peaks on detector-centroided MS1 frames have no
    /// same-m/z partner in the previous IM scan within 100 ppm). Raise it when
    /// multi-scan hills are typical, as in Biosaur2's original RT-axis use.
    ///
    /// 0 is treated as 1, as in the source: a hill has at least one peak, and
    /// the valley search would otherwise start at index -1.
    pub min_hill_length: usize,
    /// Input peaks with an intensity below this are discarded before linking.
    /// Default 0.
    pub min_intensity: f64,
    /// Tolerance for treating two IM values as the same IM scan when no scan
    /// indices are given. 0 (the default) means exact equality, which usually
    /// works because TIMS calibration yields one discrete IM value per scan;
    /// callers whose IM values went through `f32` may want a small epsilon.
    pub im_group_tolerance: f64,
    /// Maximum number of consecutive empty IM scans a hill may bridge while
    /// linking. 0 (the default) is strict consecutive-scan linking, Biosaur2's
    /// behaviour; 1 tolerates one empty scan, so the hill can still be
    /// extended two scans later. The hill length counts only the scans where
    /// the ion was observed, not the bridged gap.
    ///
    /// A hill of length 1 never bridges a gap: its tail must be extended in the
    /// very next scan or it is dropped, so two unrelated single-scan blips
    /// cannot chain into a false length-2 hill.
    ///
    /// Requires scan indices: without them empty scans cannot be told apart
    /// from absent ones, and [`centroid_frame`] uses 0, as the source does.
    pub max_scan_gap: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            mz_tol_ppm: 8.0,
            valley_factor: 1.3,
            min_hill_length: 1,
            min_intensity: 0.0,
            im_group_tolerance: 0.0,
            max_scan_gap: 0,
        }
    }
}

/// One centroided peak produced by hill detection.
///
/// The source struct `PASEFHillCentroider::Centroid`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Centroid {
    /// Intensity-weighted mean m/z over the hill.
    pub mz: f64,
    /// Summed intensity across the hill.
    pub intensity: f64,
    /// IM value at the hill's intensity apex (the first peak of maximal
    /// intensity; 0-intensity hills use their first peak).
    pub im_apex: f64,
    /// Smallest IM value covered by the hill.
    pub im_lower: f64,
    /// Largest IM value covered by the hill.
    pub im_upper: f64,
    /// Smallest m/z value covered by the hill.
    pub mz_lower: f64,
    /// Largest m/z value covered by the hill.
    pub mz_upper: f64,
    /// Number of IM scans contributing to the hill (bridged gaps excluded).
    pub length: usize,
}

/// A trace of one ion across the IM scans of a frame, arrays in scan order.
#[derive(Clone, Debug, Default)]
struct Hill {
    mz_values: Vec<f64>,
    intensities: Vec<f64>,
    im_values: Vec<f64>,
    intensity_sum: f64,
    intensity_apex: f64,
    apex_idx: usize,
}

impl Hill {
    fn push(&mut self, mz: f64, intensity: f64, im: f64) {
        self.mz_values.push(mz);
        self.intensities.push(intensity);
        self.im_values.push(im);
        self.intensity_sum += intensity;
        if intensity > self.intensity_apex {
            self.intensity_apex = intensity;
            self.apex_idx = self.intensities.len() - 1;
        }
    }
}

/// The most recent peak of an active hill and the scans since its extension.
#[derive(Clone, Copy, Debug)]
struct Tail {
    hill_id: usize,
    mz: f64,
    intensity: f64,
    fm: i64,
    /// 0: extended in the most recent scan; `n`: `n` scans without extension.
    age: usize,
}

/// Linking state carried from scan to scan.
#[derive(Debug, Default)]
struct Linker {
    hills: Vec<Hill>,
    tails: Vec<Tail>,
    /// Tail indices by m/z bin, in tail order.
    tails_by_fm: BTreeMap<i64, Vec<usize>>,
}

/// What one current-scan peak did during linking.
#[derive(Clone, Copy)]
enum Assignment {
    /// Skipped (m/z not positive).
    Skipped,
    /// Extended the tail with this index.
    Extended(usize),
    /// Started the hill with this id.
    Started(usize),
}

impl Linker {
    /// Every tail gets `scans` scans older, then the source's prune runs once.
    fn age_empty_scans(&mut self, scans: usize, max_scan_gap: usize) {
        if scans == 0 || self.tails.is_empty() {
            return;
        }
        for tail in &mut self.tails {
            tail.age = tail.age.saturating_add(scans);
        }
        self.prune_and_index(max_scan_gap);
    }

    /// The source `linkScan` for a non-empty scan, then `pruneAndIndexTails`.
    fn link_scan(
        &mut self,
        peaks: &[(f64, f64, i64)],
        scan_im: f64,
        mz_tol_ppm: f64,
        max_scan_gap: usize,
    ) {
        // Greedy assignment in descending current-scan intensity.
        let mut order: Vec<usize> = (0..peaks.len()).collect();
        order.sort_by(|&a, &b| {
            peaks[b]
                .1
                .partial_cmp(&peaks[a].1)
                .unwrap_or(Ordering::Equal)
        });
        let mut claimed = vec![false; self.tails.len()];
        let mut assignment = vec![Assignment::Skipped; peaks.len()];
        for &i in &order {
            let (mz, intensity, fm) = peaks[i];
            if mz <= 0.0 {
                continue;
            }
            let mut best_intensity = 0.0;
            let mut best_tail = None;
            for bin in [fm, fm.saturating_sub(1), fm.saturating_add(1)] {
                let Some(candidates) = self.tails_by_fm.get(&bin) else {
                    continue;
                };
                for &index in candidates {
                    if claimed[index] {
                        continue;
                    }
                    let tail = &self.tails[index];
                    let dppm = (mz - tail.mz).abs() / mz * 1e6;
                    if dppm > mz_tol_ppm {
                        continue;
                    }
                    // Biosaur2 tie-break: `>=` prefers the later tail.
                    if tail.intensity >= best_intensity {
                        best_intensity = tail.intensity;
                        best_tail = Some(index);
                    }
                }
            }
            if let Some(index) = best_tail {
                claimed[index] = true;
                assignment[i] = Assignment::Extended(index);
                let hill_id = self.tails[index].hill_id;
                self.hills[hill_id].push(mz, intensity, scan_im);
            } else {
                let hill_id = self.hills.len();
                let mut hill = Hill::default();
                hill.push(mz, intensity, scan_im);
                self.hills.push(hill);
                assignment[i] = Assignment::Started(hill_id);
            }
        }

        for tail in &mut self.tails {
            tail.age = tail.age.saturating_add(1);
        }
        for (i, &(mz, intensity, fm)) in peaks.iter().enumerate() {
            match assignment[i] {
                Assignment::Skipped => {}
                Assignment::Extended(index) => {
                    let tail = &mut self.tails[index];
                    tail.mz = mz;
                    tail.intensity = intensity;
                    tail.fm = fm;
                    tail.age = 0;
                }
                Assignment::Started(hill_id) => self.tails.push(Tail {
                    hill_id,
                    mz,
                    intensity,
                    fm,
                    age: 0,
                }),
            }
        }
        self.prune_and_index(max_scan_gap);
    }

    /// The source `pruneAndIndexTails`: drop tails older than `max_scan_gap`
    /// and length-1 tails that were not extended, then rebuild the bin index.
    fn prune_and_index(&mut self, max_scan_gap: usize) {
        let hills = &self.hills;
        self.tails.retain(|tail| {
            let length = hills[tail.hill_id].intensities.len();
            !(length < 2 && tail.age > 0) && tail.age <= max_scan_gap
        });
        self.tails_by_fm.clear();
        for (index, tail) in self.tails.iter().enumerate() {
            self.tails_by_fm.entry(tail.fm).or_default().push(index);
        }
    }
}

/// 3-tap symmetric mean filter with zero padding (Biosaur2 `meanFilter_`,
/// window 1).
fn mean_filter3(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    (0..n)
        .map(|i| {
            let mut sum = data[i];
            if i > 0 {
                sum += data[i - 1];
            }
            if i + 1 < n {
                sum += data[i + 1];
            }
            sum / 3.0
        })
        .collect()
}

/// The source `splitHills`, a port of Biosaur2's `splitHills_`.
fn split_hills(hills: Vec<Hill>, hvf: f64, min_length: usize) -> Vec<Hill> {
    let mut out = Vec::with_capacity(hills.len());
    for hill in hills {
        let hill_len = hill.intensities.len();
        // `hill_len < 2 * min_length`, without the source's `Size` wrap.
        if min_length > hill_len / 2 {
            out.push(hill);
            continue;
        }
        let smoothed = mean_filter3(&hill.intensities);
        let mut min_idx_list: Vec<usize> = Vec::new();
        let mut recheck_positions: Vec<usize> = Vec::new();
        // hill_len >= 2 * min_length >= 2, so neither underflows.
        let c_len = hill_len - min_length;
        let mut idx = min_length - 1;
        let mut l_idx = 0usize;
        let mut min_val = 0.0;
        while idx <= c_len {
            if let Some(&last) = min_idx_list.last() {
                if idx >= last + min_length {
                    l_idx = last;
                }
            }
            let valley = smoothed[idx];
            if valley <= 0.0 {
                idx += 1;
                continue;
            }
            let left_max = smoothed[l_idx..idx].iter().fold(0.0f64, |m, &v| m.max(v));
            if left_max == 0.0 {
                idx += 1;
                continue;
            }
            let l_r = left_max / valley;
            if l_r >= hvf {
                let right_max = smoothed[idx + 1..].iter().fold(0.0f64, |m, &v| m.max(v));
                if right_max > 0.0 {
                    let r_r = right_max / valley;
                    if r_r >= hvf {
                        let mult_val = l_r * r_r;
                        let candidate = idx + usize::from(l_r > r_r);
                        if min_length <= candidate && candidate <= c_len {
                            match min_idx_list.last() {
                                Some(&last) if candidate < last + min_length => {
                                    if mult_val > min_val {
                                        *min_idx_list.last_mut().expect("non-empty") = candidate;
                                        *recheck_positions.last_mut().expect("non-empty") = idx;
                                        min_val = mult_val;
                                    }
                                }
                                _ => {
                                    min_idx_list.push(candidate);
                                    recheck_positions.push(idx);
                                    min_val = mult_val;
                                }
                            }
                        }
                    }
                }
            }
            idx += 1;
        }

        // Second pass: confirm the right-side maximum for each candidate.
        let mut final_splits = Vec::new();
        for (k, &split) in min_idx_list.iter().enumerate() {
            let recheck = recheck_positions[k];
            let end = min_idx_list.get(k + 1).copied().unwrap_or(hill_len);
            if recheck + 1 >= end || recheck >= hill_len {
                continue;
            }
            let recheck_val = smoothed[recheck];
            if recheck_val <= 0.0 {
                continue;
            }
            let right_max = smoothed[recheck + 1..end]
                .iter()
                .fold(0.0f64, |m, &v| m.max(v));
            if right_max / recheck_val >= hvf {
                final_splits.push(split);
            }
        }
        if final_splits.is_empty() {
            out.push(hill);
            continue;
        }

        let mut boundaries = Vec::with_capacity(final_splits.len() + 2);
        boundaries.push(0);
        boundaries.extend_from_slice(&final_splits);
        boundaries.push(hill_len);
        for pair in boundaries.windows(2) {
            let (start, end) = (pair[0], pair[1]);
            if end <= start || end - start < min_length {
                continue;
            }
            let mut segment = Hill {
                mz_values: hill.mz_values[start..end].to_vec(),
                intensities: hill.intensities[start..end].to_vec(),
                im_values: hill.im_values[start..end].to_vec(),
                ..Hill::default()
            };
            let mut sum = 0.0;
            let mut apex = 0.0;
            let mut apex_idx = 0;
            for (k, &intensity) in segment.intensities.iter().enumerate() {
                sum += intensity;
                if intensity > apex {
                    apex = intensity;
                    apex_idx = k;
                }
            }
            segment.intensity_sum = sum;
            segment.intensity_apex = apex;
            segment.apex_idx = apex_idx;
            out.push(segment);
        }
    }
    out
}

fn require_finite(values: &[f64], what: &str) -> Result<()> {
    if let Some(index) = values.iter().position(|v| !v.is_finite()) {
        return Err(Error::InvalidValue(format!(
            "PASEF hill centroiding: {what} value {} at index {index} is not finite",
            values[index]
        )));
    }
    Ok(())
}

/// Centroid one TIMS-PASEF frame by IM-axis hill detection.
///
/// The source `centroidFrame(mz_values, intensities, im_values, n, p,
/// scan_indices)`. `mz_values`, `intensities` and `im_values` are the flat
/// peaks of one frame and need not be sorted; the function sorts internally.
///
/// When `scan_indices` is given, it holds the integer TIMS IM-scan index of
/// each peak; the linker then treats every index from the smallest to the
/// largest as a scan, so empty scans age active hill tails and a tail is closed
/// once its age exceeds [`Params::max_scan_gap`]. Without it, peaks are grouped
/// by IM value within [`Params::im_group_tolerance`] (all peaks of one scan
/// must share their IM value to within it), only non-empty groups are visited,
/// and gap bridging is disabled.
///
/// Returns the centroids sorted by ascending m/z; empty when there are no
/// peaks or every peak fails the filters.
///
/// # Errors
///
/// [`Error::InvalidValue`] when the slices differ in length, when any m/z,
/// intensity or IM value is not finite, or when [`Params::mz_tol_ppm`],
/// [`Params::valley_factor`], [`Params::min_intensity`] or
/// [`Params::im_group_tolerance`] is not finite. The source has no error path;
/// these inputs are undefined behaviour there (see the module documentation).
pub fn centroid_frame(
    mz_values: &[f64],
    intensities: &[f64],
    im_values: &[f64],
    params: &Params,
    scan_indices: Option<&[u32]>,
) -> Result<Vec<Centroid>> {
    let n = mz_values.len();
    if intensities.len() != n
        || im_values.len() != n
        || scan_indices.is_some_and(|scans| scans.len() != n)
    {
        return Err(Error::InvalidValue(
            "PASEF hill centroiding: m/z, intensity, IM and scan-index arrays differ in length"
                .into(),
        ));
    }
    for (value, name) in [
        (params.mz_tol_ppm, "mz_tol_ppm"),
        (params.valley_factor, "valley_factor"),
        (params.min_intensity, "min_intensity"),
        (params.im_group_tolerance, "im_group_tolerance"),
    ] {
        if !value.is_finite() {
            return Err(Error::InvalidValue(format!(
                "PASEF hill centroiding: parameter {name} is not finite"
            )));
        }
    }
    require_finite(mz_values, "m/z")?;
    require_finite(intensities, "intensity")?;
    require_finite(im_values, "IM")?;
    if n == 0 {
        return Ok(Vec::new());
    }
    let min_hill_length = params.min_hill_length.max(1);

    let mut max_mz = 0.0f64;
    let mut kept = Vec::with_capacity(n);
    for i in 0..n {
        if intensities[i] < params.min_intensity {
            continue;
        }
        if mz_values[i] > max_mz {
            max_mz = mz_values[i];
        }
        kept.push(i);
    }
    if kept.is_empty() {
        return Ok(Vec::new());
    }
    let mz_step = (params.mz_tol_ppm * 1e-6) * max_mz;
    if mz_step <= 0.0 {
        return Ok(Vec::new());
    }
    // The source's `static_cast<int>(mz / mz_step)`: truncation toward zero.
    // The cast saturates where the source is undefined.
    let bin = |mz: f64| (mz / mz_step) as i64;
    let by_mz = |a: usize, b: usize| {
        mz_values[a]
            .partial_cmp(&mz_values[b])
            .unwrap_or(Ordering::Equal)
    };

    let mut linker = Linker::default();
    let mut peaks: Vec<(f64, f64, i64)> = Vec::new();
    if let Some(scans) = scan_indices {
        kept.sort_by(|&a, &b| scans[a].cmp(&scans[b]).then_with(|| by_mz(a, b)));
        let mut previous: Option<u32> = None;
        let mut start = 0;
        while start < kept.len() {
            let scan = scans[kept[start]];
            let end = start
                + kept[start..]
                    .iter()
                    .take_while(|&&i| scans[i] == scan)
                    .count();
            if let Some(previous) = previous {
                let empty = (scan - previous - 1) as usize;
                linker.age_empty_scans(empty, params.max_scan_gap);
            }
            peaks.clear();
            peaks.extend(
                kept[start..end]
                    .iter()
                    .map(|&i| (mz_values[i], intensities[i], bin(mz_values[i]))),
            );
            // All peaks of one scan share their IM value; the source keeps
            // the last one's.
            let scan_im = im_values[kept[end - 1]];
            linker.link_scan(&peaks, scan_im, params.mz_tol_ppm, params.max_scan_gap);
            previous = Some(scan);
            start = end;
        }
    } else {
        kept.sort_by(|&a, &b| {
            im_values[a]
                .partial_cmp(&im_values[b])
                .unwrap_or(Ordering::Equal)
                .then_with(|| by_mz(a, b))
        });
        let mut start = 0;
        while start < kept.len() {
            let scan_im = im_values[kept[start]];
            let end = start
                + kept[start..]
                    .iter()
                    .take_while(|&&i| (im_values[i] - scan_im).abs() <= params.im_group_tolerance)
                    .count();
            peaks.clear();
            peaks.extend(
                kept[start..end]
                    .iter()
                    .map(|&i| (mz_values[i], intensities[i], bin(mz_values[i]))),
            );
            // Gap bridging is forced off without scan indices.
            linker.link_scan(&peaks, scan_im, params.mz_tol_ppm, 0);
            start = end;
        }
    }

    let hills = split_hills(linker.hills, params.valley_factor, min_hill_length);
    let mut result = Vec::with_capacity(hills.len());
    for hill in hills {
        let length = hill.intensities.len();
        if length < min_hill_length {
            continue;
        }
        let mut weight = 0.0;
        let mut weighted_mz = 0.0;
        let mut im_lower = hill.im_values[0];
        let mut im_upper = hill.im_values[0];
        let mut mz_lower = hill.mz_values[0];
        let mut mz_upper = hill.mz_values[0];
        for k in 0..length {
            let (mz, intensity, im) = (hill.mz_values[k], hill.intensities[k], hill.im_values[k]);
            weight += intensity;
            weighted_mz += mz * intensity;
            if im < im_lower {
                im_lower = im;
            }
            if im > im_upper {
                im_upper = im;
            }
            if mz < mz_lower {
                mz_lower = mz;
            }
            if mz > mz_upper {
                mz_upper = mz;
            }
        }
        if weight <= 0.0 {
            continue;
        }
        result.push(Centroid {
            mz: weighted_mz / weight,
            intensity: hill.intensity_sum,
            im_apex: hill.im_values[hill.apex_idx],
            im_lower,
            im_upper,
            mz_lower,
            mz_upper,
            length,
        });
    }
    result.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_filter_pads_with_zero() {
        assert_eq!(mean_filter3(&[3.0, 6.0, 9.0]), vec![3.0, 6.0, 5.0]);
        assert_eq!(mean_filter3(&[3.0]), vec![1.0]);
    }

    #[test]
    fn mismatched_lengths_and_nan_are_refused() {
        let params = Params::default();
        assert!(centroid_frame(&[1.0], &[], &[1.0], &params, None).is_err());
        assert!(centroid_frame(&[1.0], &[1.0], &[1.0], &params, Some(&[])).is_err());
        assert!(centroid_frame(&[f64::NAN], &[1.0], &[1.0], &params, None).is_err());
        let params = Params {
            mz_tol_ppm: f64::NAN,
            ..Params::default()
        };
        assert!(centroid_frame(&[1.0], &[1.0], &[1.0], &params, None).is_err());
    }

    #[test]
    fn largest_scan_index_terminates() {
        // The source loop `for (s = scan_min; s <= scan_max; ++s)` never ends
        // for scan_max == u32::MAX.
        let result = centroid_frame(
            &[500.0, 500.0],
            &[100.0, 200.0],
            &[1.0, 1.1],
            &Params::default(),
            Some(&[u32::MAX - 1, u32::MAX]),
        )
        .unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].length, 2);
    }

    #[test]
    fn zero_min_hill_length_is_one() {
        let params = Params {
            min_hill_length: 0,
            ..Params::default()
        };
        let result = centroid_frame(&[500.0], &[10.0], &[1.0], &params, None).unwrap();
        assert_eq!(result.len(), 1);
    }
}
