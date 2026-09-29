// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The remainder of `CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h` on
//! [`IsotopeDistribution`]: `merge`, `operator<`, `std::hash`, the `Sorted`
//! enum and validated in-place access standing in for the mutable iterators
//! and `operator[]`.
//!
//! # API mapping
//!
//! | Source (`IsotopeDistribution.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `enum Sorted {INTENSITY, MASS, UNDEFINED}` | [`Sorted`] |
//! | `merge(double resolution, double min_prob)` | [`IsotopeDistribution::merge`] |
//! | `operator<` | `Ord` |
//! | `std::hash<IsotopeDistribution>` | `Hash`; the source's value is [`IsotopeDistribution::source_hash64`] |
//! | `Peak1D& operator[](Size)` | [`IsotopeDistribution::set_peak`] (validated) |
//! | `Iterator begin()`/`end()`, `ReverseIterator rbegin()`/`rend()` (mutable) | [`IsotopeDistribution::try_for_each_peak_mut`] (validated) |
//! | protected `sort_`, `transform_` | the sorts and [`IsotopeDistribution::try_for_each_peak_mut`] |

use super::{IsotopeDistribution, IsotopePeak, invalid, validate_nonnegative};
use crate::Result;
use crate::concept::hash_utils::{hash_combine, hash_float, hash_int};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

/// Source `IsotopeDistribution::Sorted`: the order a distribution is in.
///
/// The source declares the enum but no member reads or stores it; it is kept
/// for callers that name the three states.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Sorted {
    /// Source `INTENSITY`: descending probability.
    Intensity,
    /// Source `MASS`: ascending mass.
    Mass,
    /// Source `UNDEFINED`.
    #[default]
    Undefined,
}

impl IsotopeDistribution {
    /// Re-bins the distribution onto a coarser mass grid: source
    /// `merge(double resolution, double min_prob)`.
    ///
    /// Following `IsotopeDistribution.cpp:240-273`: sort by mass, trim both
    /// tails below `min_prob` (the left trim keeps every peak when none reaches
    /// it, as source `trimLeft` does), then allocate
    /// `n = ceil((max - min) / resolution)` bins spaced `(max - min) / n` from
    /// the lightest mass. Each peak goes to bin
    /// `round((mass - min) / resolution)`; a peak whose bin is `n` or beyond is
    /// dropped, and a bin's mass is overwritten by every peak landing in it
    /// while its weights add. Finally weights below `min_prob` are removed,
    /// including bins nothing landed in, which keep mass and weight zero when
    /// `min_prob` is zero.
    ///
    /// Two consequences of that rule are the source's and kept: when
    /// `max - min` is an exact multiple of `resolution` the heaviest peak's bin
    /// is `n` and its weight is dropped, and a distribution whose remaining
    /// peaks share one mass (`n = 0`) becomes empty.
    ///
    /// Weights add in `f64`; the source adds them in the `float` intensity of
    /// `Peak1D`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when `resolution` is not positive
    /// and finite, `min_prob` is negative or not finite, the trimmed
    /// distribution is empty (the source reads past its end), or `n` exceeds
    /// the number of remaining peaks (the source's `IllegalArgument`, "New
    /// Isotope Distribution has more points than the old one."). The
    /// distribution is then unchanged.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::isotopes::{IsotopeDistribution, IsotopePeak};
    ///
    /// let peak = |mass, probability| IsotopePeak { mass, probability };
    /// let mut distribution = IsotopeDistribution::from_peaks(vec![
    ///     peak(100.0, 0.5),
    ///     peak(100.4, 0.25),
    ///     peak(101.0, 0.125),
    ///     peak(101.5, 0.125),
    /// ])
    /// .unwrap();
    /// distribution.merge(1.0, 0.0).unwrap();
    /// // Two bins of width 0.75 from 100.0; the peak at 101.5 rounds to bin 2
    /// // and is dropped.
    /// let merged: Vec<_> = distribution.peaks().iter().map(|p| (p.mass, p.probability)).collect();
    /// assert_eq!(merged, [(100.0, 0.75), (100.75, 0.125)]);
    /// ```
    pub fn merge(&mut self, resolution: f64, min_prob: f64) -> Result<()> {
        if !resolution.is_finite() || resolution <= 0.0 {
            return Err(invalid("merge resolution must be finite and positive"));
        }
        validate_nonnegative(min_prob, "merge probability cutoff")?;
        let mut raw = self.clone();
        raw.sort_by_mass();
        // Source trimLeft: erase before the first weight reaching the cutoff,
        // nothing when none does.
        if let Some(first) = raw.peaks.iter().position(|p| p.probability >= min_prob) {
            raw.peaks.drain(..first);
        }
        raw.trim_right(min_prob)?;
        let (Some(front), Some(back)) = (raw.peaks.first(), raw.peaks.last()) else {
            return Err(invalid(
                "merge needs at least one peak at or above the probability cutoff",
            ));
        };
        let front = front.mass;
        let range = back.mass - front;
        let bins = (range / resolution).ceil();
        if bins > raw.peaks.len() as f64 {
            return Err(invalid(
                "New Isotope Distribution has more points than the old one.",
            ));
        }
        // bins <= raw.len() <= MAX_ISOTOPE_PEAKS, so the conversion is exact.
        let bins = bins as usize;
        let delta = range / bins as f64;
        let mut merged = vec![IsotopePeak::default(); bins];
        for peak in &raw.peaks {
            let index = ((peak.mass - front) / resolution).round();
            if index >= bins as f64 {
                continue;
            }
            let index = index as usize;
            merged[index].mass = front + index as f64 * delta;
            merged[index].probability += peak.probability;
        }
        for peak in &merged {
            if !peak.probability.is_finite() {
                return Err(invalid("merged isotope probability overflows f64"));
            }
        }
        merged.retain(|p| p.probability >= min_prob);
        self.peaks = merged;
        Ok(())
    }

    /// Replaces the peak at `index`: the validated form of source
    /// `Peak1D& operator[](Size)`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when `index` is out of range or
    /// the mass or probability is negative or not finite; the source indexes
    /// without a check. The distribution is then unchanged.
    pub fn set_peak(&mut self, index: usize, peak: IsotopePeak) -> Result<()> {
        validate_nonnegative(peak.mass, "isotope mass")?;
        validate_nonnegative(peak.probability, "isotope probability")?;
        let slot = self
            .peaks
            .get_mut(index)
            .ok_or_else(|| invalid(format!("isotope peak index {index} out of range")))?;
        *slot = peak;
        Ok(())
    }

    /// Applies `edit` to every peak in container order and keeps the result
    /// only when every edited peak is still valid: the validated form of the
    /// source's mutable `begin()`/`end()` and `rbegin()`/`rend()` iteration and
    /// its protected `transform_`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when an edited mass or probability
    /// is negative or not finite; the distribution is then unchanged.
    pub fn try_for_each_peak_mut<F>(&mut self, mut edit: F) -> Result<()>
    where
        F: FnMut(&mut IsotopePeak),
    {
        let mut peaks = self.peaks.clone();
        for peak in &mut peaks {
            edit(peak);
            validate_nonnegative(peak.mass, "isotope mass")?;
            validate_nonnegative(peak.probability, "isotope probability")?;
        }
        self.peaks = peaks;
        Ok(())
    }

    /// The value source `std::hash<IsotopeDistribution>` computes.
    ///
    /// The seed is `hash_int` of the peak count as the eight bytes of `Size`;
    /// each peak's `std::hash<Peak1D>` (`hash_float` of the `double` mass,
    /// combined with `hash_float` of the `float` intensity) is combined into it.
    /// The probability is narrowed to `f32` first, as the source's `Peak1D`
    /// stores it; a weight beyond the binary32 range narrows to infinity, as a
    /// source `float` assignment does.
    pub fn source_hash64(&self) -> u64 {
        let mut seed = hash_int(self.peaks.len() as u64);
        for peak in &self.peaks {
            let mut peak_hash = hash_float(peak.mass);
            hash_combine(&mut peak_hash, hash_float(peak.probability as f32));
            hash_combine(&mut seed, peak_hash);
        }
        seed
    }
}

impl Eq for IsotopeDistribution {}

impl PartialOrd for IsotopeDistribution {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IsotopeDistribution {
    /// Source `operator<`: fewer peaks first; at equal size the first differing
    /// peak decides by mass, then by probability.
    ///
    /// Masses and probabilities are finite, so every pair compares; `-0.0` and
    /// `0.0` are equal, as they are for `==`.
    fn cmp(&self, other: &Self) -> Ordering {
        self.peaks.len().cmp(&other.peaks.len()).then_with(|| {
            for (left, right) in self.peaks.iter().zip(&other.peaks) {
                let order = compare(left.mass, right.mass)
                    .then_with(|| compare(left.probability, right.probability));
                if order != Ordering::Equal {
                    return order;
                }
            }
            Ordering::Equal
        })
    }
}

impl Hash for IsotopeDistribution {
    /// Consistent with `==`: `-0.0` hashes as `0.0`.
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.peaks.len().hash(state);
        for peak in &self.peaks {
            normalized_bits(peak.mass).hash(state);
            normalized_bits(peak.probability).hash(state);
        }
    }
}

fn compare(left: f64, right: f64) -> Ordering {
    left.partial_cmp(&right).unwrap_or(Ordering::Equal)
}

fn normalized_bits(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}
