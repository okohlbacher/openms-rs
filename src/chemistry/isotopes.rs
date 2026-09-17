// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Coarse isotope distributions and averagine estimates.
//!
//! Ported from `CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h`,
//! `CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h`, their
//! implementations and the averagine estimators of `CHEMISTRY/EmpiricalFormula.cpp`.
//! The native port was reviewed at OpenMS4-core `7c029e8`; both implementation
//! files and both headers are hash-identical at the current pin `bc9cc12`.
//!
//! Probabilities are stored as `f64`. By default ([`ProbabilityPrecision::Double`])
//! they are also computed in `f64`, whereas C++ keeps them in `Peak1D`'s `float`
//! intensity. [`ProbabilityPrecision::SourceSingle`] selects the source binary32
//! arithmetic explicitly. For formulas of natural elements it reproduces, bit for
//! bit, the executed C++ SDK runs whose element iteration order is ascending
//! atomic number, which were the majority of the measured runs, except formulas
//! containing iridium, which the SDK's `ElementDB` builds from rhenium's tables
//! (`ElementDB.cpp:512`). The SDK's own order, and with it its binary32 output,
//! varies between runs. Coarse mass correction is a carbon-13 spacing
//! approximation, not isotope fine structure.
//! See `docs/ISOTOPE_SUPPORT.md` for the native conventions and limits and
//! `docs/ISOTOPE_SOURCE_PRECISION_SUPPORT.md` for the source-precision contract.

use super::{Atom, C13C12_MASSDIFF_U, EmpiricalFormula, PROTON_MASS_U, element, parse_atom};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Maximum length of any allocated isotope vector, including gap filling.
pub const MAX_ISOTOPE_PEAKS: usize = 1_000_000;
/// Maximum probability products in one coarse run or raw convolution/power call.
pub const MAX_CONVOLUTION_PRODUCTS: usize = 50_000_000;
const MAX_EXACT_INTEGER: f64 = 9_007_199_254_740_991.0;
const NEUTRON_MASS_U: f64 = 1.008_664_915_66;

/// A mass in daltons and its nonnegative probability or unnormalized weight.
///
/// Source `IsotopeDistribution::MassAbundance` is a `Peak1D`: the mass sits in
/// its `double` m/z slot and the probability in its `float` intensity slot. Both
/// are `f64` here; values produced under [`ProbabilityPrecision::SourceSingle`]
/// are exactly representable as `f32`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsotopePeak {
    /// Isotope or coarse-bin mass in daltons (the source m/z slot).
    pub mass: f64,
    /// Probability or unnormalized weight (the source intensity slot).
    pub probability: f64,
}

/// Floating-point contract for coarse isotope probabilities.
///
/// C++ stores isotope probabilities in `Peak1D`'s `float` intensity, so every
/// convolution step and the final renormalization of the source run in
/// binary32. The native default computes in `f64`; it differs from the source
/// by rounding (about 1e-7 relative). Code that must reproduce the source's
/// exact values, for example thresholds applied to binary32 patterns, selects
/// [`Self::SourceSingle`] explicitly. Narrowing is never the default because it
/// discards precision.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ProbabilityPrecision {
    /// Native `f64` arithmetic throughout; the unchanged default.
    #[default]
    Double,
    /// Source `Peak1D` binary32 arithmetic, as executed by the C++ SDK.
    ///
    /// - Input weights (element tables and isotope overrides) are narrowed to
    ///   `f32` when read, as a `float` intensity assignment does: round to
    ///   nearest, ties to even. A weight beyond the binary32 range is an error.
    /// - Each convolution step evaluates `p + l * r` in `f32` as two separate
    ///   operations. A C++ build that contracts the expression into a fused
    ///   multiply-add would differ in the last bit; the executed SDK build does
    ///   not.
    /// - Formula elements are convolved in ascending atomic number, each natural
    ///   element before its labelled isotopes in ascending mass number, and the
    ///   lightest-isotope mass is summed in the same order. The source iterates
    ///   `EmpiricalFormula`'s `std::map<const Element*, SignedSize>`
    ///   (`EmpiricalFormula.h:66`), so its order is the address order of the
    ///   `ElementDB` elements, and that order is not fixed. In 200 runs of one
    ///   SDK binary, `C1H1N1O1S1P1` iterated `H C N O P S` in 198 runs and
    ///   `H N C O P S` in 2; those 2 runs produced different binary32 patterns,
    ///   including every FeatureFinderAlgorithmPicked averagine window from
    ///   150 to 8050 Da. No fixed order reproduces every run. Ascending atomic
    ///   number reproduces the majority runs for natural elements, except
    ///   formulas containing iridium, which the SDK's `ElementDB` builds from
    ///   rhenium's tables (`ElementDB.cpp:512`). Labelled isotopes had no
    ///   majority position, so their placement here is a native choice, and the
    ///   lightest-isotope mass of a labelled formula can differ from a C++ run in
    ///   the last bit.
    /// - Renormalization sums the binary32 weights in reverse order into an
    ///   `f64` and narrows each quotient back to `f32`. When every retained bin
    ///   has underflowed to zero, as for a 1,000,000 Da peptide averagine
    ///   estimate limited to 20 bins, this returns an error where
    ///   [`Self::Double`] succeeds; the source divides zero by zero and returns
    ///   NaN weights.
    ///
    /// Results are stored as `f64` values that are exactly representable in
    /// `f32`.
    SourceSingle,
}

/// Validated isotope masses and weights, preserving insertion order.
///
/// Source `IsotopeDistribution` is a container of masses and their
/// probabilities; the calculations are done by pattern generators such as
/// [`CoarseIsotopePatternGenerator`] or the fine isotope generator. The
/// default is the convolution identity `(0, 1)`, as in the source default
/// constructor. Use [`Self::empty`] for an empty container. Mutating methods
/// preserve finite, nonnegative weights and masses; read-only peak slices cannot
/// bypass validation.
#[derive(Clone, Debug, PartialEq)]
pub struct IsotopeDistribution {
    peaks: Vec<IsotopePeak>,
}

impl Default for IsotopeDistribution {
    fn default() -> Self {
        Self {
            peaks: vec![IsotopePeak {
                mass: 0.0,
                probability: 1.0,
            }],
        }
    }
}

impl IsotopeDistribution {
    /// An empty distribution with no peaks.
    ///
    /// The source has no empty constructor; clearing a default-constructed
    /// distribution gives the same state.
    pub fn empty() -> Self {
        Self { peaks: Vec::new() }
    }

    /// A distribution holding `peaks` in the given order (source `set`).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when there are more than
    /// [`MAX_ISOTOPE_PEAKS`] peaks or a mass or probability is negative or not
    /// finite. The source `set` replaces the container without checks.
    pub fn from_peaks(peaks: Vec<IsotopePeak>) -> Result<Self> {
        check_size(peaks.len())?;
        for peak in &peaks {
            validate_nonnegative(peak.mass, "isotope mass")?;
            validate_nonnegative(peak.probability, "isotope probability")?;
        }
        Ok(Self { peaks })
    }

    /// The peaks in container order (source `getContainer`, iterators and
    /// `operator[]`).
    pub fn peaks(&self) -> &[IsotopePeak] {
        &self.peaks
    }
    /// The number of peaks (source `size`).
    pub fn len(&self) -> usize {
        self.peaks.len()
    }
    /// Whether the distribution holds no peaks.
    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty()
    }
    /// Remove every peak (source `clear`).
    pub fn clear(&mut self) {
        self.peaks.clear();
    }

    /// Append a peak without sorting (source `insert(mass, intensity)`).
    ///
    /// The source narrows the intensity to `float` on insertion; this stores the
    /// `f64` value, and [`ProbabilityPrecision::SourceSingle`] operations narrow
    /// it when they read it, which yields the same binary32 value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the distribution already holds
    /// [`MAX_ISOTOPE_PEAKS`] peaks or the mass or probability is negative or not
    /// finite; the distribution is then unchanged.
    pub fn insert(&mut self, peak: IsotopePeak) -> Result<()> {
        check_size(self.len().saturating_add(1))?;
        validate_nonnegative(peak.mass, "isotope mass")?;
        validate_nonnegative(peak.probability, "isotope probability")?;
        self.peaks.push(peak);
        Ok(())
    }

    /// Resize, padding with `(0, 0)` as in OpenMS. Does not preserve sortedness.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `length` exceeds
    /// [`MAX_ISOTOPE_PEAKS`]; the distribution is then unchanged.
    pub fn resize(&mut self, length: usize) -> Result<()> {
        check_size(length)?;
        self.peaks.resize(length, IsotopePeak::default());
        Ok(())
    }

    /// The smallest mass (source `getMin`).
    ///
    /// The source returns `0` for an empty distribution; this returns `None`.
    pub fn min_mass(&self) -> Option<f64> {
        self.peaks.iter().map(|p| p.mass).min_by(f64::total_cmp)
    }

    /// The largest mass (source `getMax`).
    ///
    /// The source returns `0` for an empty distribution; this returns `None`.
    pub fn max_mass(&self) -> Option<f64> {
        self.peaks.iter().map(|p| p.mass).max_by(f64::total_cmp)
    }

    /// Most abundant peak; the first wins ties. Empty distributions return None.
    ///
    /// Source `getMostAbundant` also keeps the first maximum
    /// (`std::max_element`), but returns `Peak1D(0, 1)` when empty.
    pub fn most_abundant(&self) -> Option<IsotopePeak> {
        self.peaks.iter().copied().reduce(|best, peak| {
            if peak.probability > best.probability {
                peak
            } else {
                best
            }
        })
    }

    /// Sort by ascending mass (source `sortByMass`).
    ///
    /// The source `std::sort` is unstable; this sort is stable, so equal masses
    /// keep their insertion order.
    pub fn sort_by_mass(&mut self) {
        self.peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
    }

    /// Sort by descending probability (source `sortByIntensity`).
    ///
    /// The source `std::sort` is unstable; this sort is stable, so equal
    /// probabilities keep their insertion order.
    pub fn sort_by_probability(&mut self) {
        self.peaks
            .sort_by(|a, b| b.probability.total_cmp(&a.probability));
    }

    /// Sum weights from the end, following the upstream accumulation direction.
    pub fn probability_sum(&self) -> f64 {
        self.peaks.iter().rev().map(|p| p.probability).sum()
    }

    /// Normalize existing weights. Empty is a no-op; zero/overflowed sums error.
    /// Failed normalization leaves the distribution unchanged.
    ///
    /// Source `renormalize` re-normalizes the sum of the probabilities of all
    /// isotopes to one; the source notes this may be needed because
    /// distributions with many isotopes accumulate inexact sums. This method
    /// computes in `f64`; [`Self::renormalize_with`] selects the source binary32
    /// contract.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the probability sum of a nonempty
    /// distribution is zero or not finite. The source divides regardless and
    /// produces NaN or infinite weights.
    pub fn renormalize(&mut self) -> Result<()> {
        if self.is_empty() {
            return Ok(());
        }
        let sum = self.probability_sum();
        if !sum.is_finite() || sum <= 0.0 {
            return Err(invalid(
                "normalization requires a finite, positive probability sum",
            ));
        }
        for peak in &mut self.peaks {
            peak.probability /= sum;
        }
        Ok(())
    }

    /// Normalize existing weights under an explicit floating-point contract.
    ///
    /// [`ProbabilityPrecision::Double`] is [`Self::renormalize`].
    /// [`ProbabilityPrecision::SourceSingle`] follows
    /// `IsotopeDistribution.cpp:176-192` on `Peak1D` intensities: every weight is
    /// narrowed to `f32`, the narrowed weights are summed from the end into an
    /// `f64`, and each weight becomes the `f32` narrowing of the `f64` quotient.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a weight exceeds the binary32 range
    /// (source-precision only) or the sum of a nonempty distribution is zero or
    /// not finite. Failed normalization leaves the distribution unchanged.
    pub fn renormalize_with(&mut self, precision: ProbabilityPrecision) -> Result<()> {
        match precision {
            ProbabilityPrecision::Double => self.renormalize(),
            ProbabilityPrecision::SourceSingle => {
                if self.is_empty() {
                    return Ok(());
                }
                let narrowed = self
                    .peaks
                    .iter()
                    .map(|peak| narrow(peak.probability))
                    .collect::<Result<Vec<f32>>>()?;
                let sum = narrowed
                    .iter()
                    .rev()
                    .fold(0.0_f64, |sum, &weight| sum + f64::from(weight));
                if !sum.is_finite() || sum <= 0.0 {
                    return Err(invalid(
                        "normalization requires a finite, positive probability sum",
                    ));
                }
                for (peak, weight) in self.peaks.iter_mut().zip(narrowed) {
                    peak.probability = f64::from((f64::from(weight) / sum) as f32);
                }
                Ok(())
            }
        }
    }

    /// Probability-weighted mass, without requiring prior normalization.
    ///
    /// Source `averageMass` divides each weight by the total before weighting,
    /// in container order; so does this.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the probability sum is zero or not
    /// finite (the source returns NaN) or the weighted mass overflows.
    pub fn average_mass(&self) -> Result<f64> {
        let sum = self.probability_sum();
        if !sum.is_finite() || sum <= 0.0 {
            return Err(invalid(
                "average mass requires a finite, positive probability sum",
            ));
        }
        let mass: f64 = self
            .peaks
            .iter()
            .map(|p| p.mass * (p.probability / sum))
            .sum();
        if !mass.is_finite() {
            return Err(invalid("weighted isotope mass overflowed"));
        }
        Ok(mass)
    }

    /// Remove all weights strictly below cutoff; equality is retained.
    ///
    /// Source `trimIntensities` removes intensities below the cutoff.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `cutoff` is negative or not finite;
    /// the distribution is then unchanged.
    pub fn trim_intensities(&mut self, cutoff: f64) -> Result<()> {
        validate_nonnegative(cutoff, "probability cutoff")?;
        self.peaks.retain(|p| p.probability >= cutoff);
        Ok(())
    }

    /// Remove the low-probability prefix in current order, without normalization.
    /// Unlike the source, this can remove every peak.
    ///
    /// Source `trimLeft` trims the left side of a distribution to isotopes with
    /// a significant contribution, typically the small leading entries of
    /// distributions calculated for large masses; the source notes that the
    /// distribution should be normalized afterwards. When no weight reaches the
    /// cutoff the source leaves the distribution unchanged; use
    /// [`Self::trim_left_source`] for that behavior.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `cutoff` is negative or not finite;
    /// the distribution is then unchanged.
    pub fn trim_left(&mut self, cutoff: f64) -> Result<()> {
        validate_nonnegative(cutoff, "probability cutoff")?;
        let first = self
            .peaks
            .iter()
            .position(|p| p.probability >= cutoff)
            .unwrap_or(self.len());
        self.peaks.drain(..first);
        Ok(())
    }

    /// Remove the low-probability prefix exactly as source `trimLeft` does.
    ///
    /// `IsotopeDistribution.cpp:210-220` erases the entries before the first
    /// weight at or above `cutoff` and erases nothing when no weight reaches it,
    /// so a distribution entirely below the cutoff is kept whole. That case is
    /// reachable in `FeatureFinderAlgorithmPicked` through
    /// `isotopic_pattern:intensity_percentage_optional`. Weights equal to the
    /// cutoff are retained. No normalization is applied.
    ///
    /// The source compares its `float` intensity, promoted to `double`, with
    /// the `double` cutoff. Each stored weight is narrowed to `f32` first, as
    /// source `insert` narrows it, so the comparison is the source's for every
    /// input: a weight of `0.7` becomes `0.699999988` and does not reach a
    /// cutoff of `0.7`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `cutoff` is negative or not finite,
    /// or a weight exceeds the binary32 range (the source stores an infinite
    /// intensity); the distribution is then unchanged.
    pub fn trim_left_source(&mut self, cutoff: f64) -> Result<()> {
        validate_nonnegative(cutoff, "probability cutoff")?;
        let narrowed = self
            .peaks
            .iter()
            .map(|peak| narrow(peak.probability))
            .collect::<Result<Vec<f32>>>()?;
        if let Some(first) = narrowed
            .iter()
            .position(|&weight| f64::from(weight) >= cutoff)
        {
            self.peaks.drain(..first);
        }
        Ok(())
    }

    /// Remove the low-probability suffix in current order, without normalization.
    ///
    /// Source `trimRight` trims the right side to isotopes with a significant
    /// contribution and notes that the distribution should be normalized
    /// afterwards. As in the source, every peak is removed when no weight
    /// reaches the cutoff, and weights equal to the cutoff are retained.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `cutoff` is negative or not finite;
    /// the distribution is then unchanged.
    pub fn trim_right(&mut self, cutoff: f64) -> Result<()> {
        validate_nonnegative(cutoff, "probability cutoff")?;
        let end = self
            .peaks
            .iter()
            .rposition(|p| p.probability >= cutoff)
            .map_or(0, |i| i + 1);
        self.peaks.truncate(end);
        Ok(())
    }
}

/// The result of [`CoarseIsotopePatternGenerator::estimate_from_peptide_weight_source`].
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum SourceSingleEstimate {
    /// The estimate, normalised in the source's binary32 arithmetic.
    Normalized(IsotopeDistribution),
    /// Every one of the `len` retained binary32 bins underflowed to zero, so
    /// the source's `renormalize` made every weight NaN (`0 / 0`).
    AllUnderflowed {
        /// The number of weights.
        len: usize,
    },
}

/// Mass labels for the same nominal-isotope probabilities.
///
/// Source `setRoundMasses(false)` (the default) corresponds to
/// [`Self::Approximate`] and `setRoundMasses(true)` to [`Self::Nominal`]. The
/// probabilities are the same coarse distribution in both modes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CoarseMassMode {
    /// Lightest-isotope mass plus `i * C13C12_MASSDIFF_U` (not exact fine mass).
    #[default]
    Approximate,
    /// Round the approximate corrected mass of each peak to the nearest integer.
    Nominal,
}

/// Nominal isotope convolution with optional low-mass-tail truncation.
///
/// Source `CoarseIsotopePatternGenerator` produces theoretical distributions for
/// empirical formulas at a resolution of 1 Da. It convolves the natural
/// abundances of each element, so the probabilities of each coarse peak are
/// accurate, but it treats every isotope as having an integer nominal mass: it
/// does not discriminate between carbon-13, nitrogen-15 and oxygen-18
/// contributions. The masses are therefore only approximately accurate; the fine
/// isotope generator resolves them. Mass labels follow [`CoarseMassMode`]; for
/// C100 the source documents
///
/// ```text
/// Nominal:      1200 : 0.341036528   1201       : 0.368855864   ...
/// Approximate:  1200 : 0.341036528   1201.00335 : 0.368855864   ...
/// ```
///
/// The maximum isotope count is an upper bound on the reported bins: three keeps
/// the monoisotopic bin, +1 and +2. The source notes that by default all
/// possible isotopes are calculated, which yields many values for large masses;
/// here that unbounded mode is additionally limited by [`MAX_ISOTOPE_PEAKS`] and
/// [`MAX_CONVOLUTION_PRODUCTS`].
///
/// `run` normalizes retained probabilities to one. `convolve`, `convolve_power`,
/// and `calc_fragment_isotope_dist` do not normalize. None means all peaks,
/// subject to documented allocation/work limits; Some(0) is invalid.
/// [`ProbabilityPrecision`] selects native `f64` or source binary32 arithmetic
/// for all of them.
#[derive(Clone, Debug, Default)]
pub struct CoarseIsotopePatternGenerator {
    max_peaks: Option<usize>,
    mass_mode: CoarseMassMode,
    precision: ProbabilityPrecision,
    overrides: BTreeMap<Atom, IsotopeDistribution>,
}

impl CoarseIsotopePatternGenerator {
    /// A generator retaining at most `max_peaks` coarse bins, in native `f64`
    /// precision.
    ///
    /// Source `CoarseIsotopePatternGenerator(max_isotope = 0, round_masses =
    /// false)`: `max_isotope` 0 means all isotopes and corresponds to `None`,
    /// and `round_masses` corresponds to `mass_mode`. Select the source binary32
    /// arithmetic with [`Self::with_precision`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for `Some(0)` (use `None` for the source's
    /// 0) and for a count above [`MAX_ISOTOPE_PEAKS`].
    pub fn new(max_peaks: Option<usize>, mass_mode: CoarseMassMode) -> Result<Self> {
        if let Some(length) = max_peaks {
            if length == 0 {
                return Err(invalid("maximum isotope count must be positive"));
            }
            check_size(length)?;
        }
        Ok(Self {
            max_peaks,
            mass_mode,
            precision: ProbabilityPrecision::Double,
            overrides: BTreeMap::new(),
        })
    }

    /// This generator with the given probability arithmetic.
    ///
    /// [`ProbabilityPrecision::SourceSingle`] reproduces the source `Peak1D`
    /// binary32 results of `run`, every `estimate_*` method, `convolve`,
    /// `convolve_power`, `calc_fragment_isotope_dist` and
    /// `estimate_fragment_from_weights`, within the element-order scope stated
    /// on that variant and except formulas containing iridium, which the SDK's
    /// `ElementDB` builds from rhenium's tables (`ElementDB.cpp:512`). The
    /// static Poisson approximations are unaffected.
    pub fn with_precision(mut self, precision: ProbabilityPrecision) -> Self {
        self.precision = precision;
        self
    }

    /// Replace the probability arithmetic; see [`Self::with_precision`].
    pub fn set_precision(&mut self, precision: ProbabilityPrecision) {
        self.precision = precision;
    }

    /// The probability arithmetic this generator uses.
    pub fn precision(&self) -> ProbabilityPrecision {
        self.precision
    }

    /// The maximum number of reported isotopes, `None` for all (source
    /// `getMaxIsotope`, where 0 means all).
    ///
    /// The source documents the limit as useful because distributions with
    /// numerous isotopes tend to end in many numerical zeros. The limit is fixed
    /// at construction; source `setMaxIsotope` corresponds to constructing a new
    /// generator.
    pub fn max_peaks(&self) -> Option<usize> {
        self.max_peaks
    }
    /// How output masses are labelled (source `getRoundMasses`).
    pub fn mass_mode(&self) -> CoarseMassMode {
        self.mass_mode
    }

    /// Override isotope weights locally, without mutating the shared element table.
    /// Masses may be exact or nominal but must round to declared isotope numbers.
    /// Include the lightest declared isotope, even at zero probability, so the
    /// mass origin remains defined. Weights need not sum to one.
    ///
    /// Source `setIsotopeOverride` uses the given distribution in place of the
    /// element's natural one wherever the element occurs in the formula, and for
    /// the implicit H+ adduct when the element is hydrogen, which keeps labelled
    /// pattern computations local and thread-safe. The source keys the override
    /// by an `ElementDB` element pointer; this takes the element or isotope
    /// symbol, and explicit isotope symbols such as `(13)C` are separate keys.
    /// Under [`ProbabilityPrecision::SourceSingle`] the weights are narrowed to
    /// `f32` when read, matching the `float` the source stores on insertion.
    ///
    /// `FeatureFinderAlgorithmPicked.cpp:166-178` builds its abundance overrides
    /// by inserting into a default-constructed source distribution, which already
    /// holds `(0, 1)`; such an input is rejected here because its origin is not
    /// the lightest declared isotope.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `symbol` is not exactly one element
    /// or isotope symbol, the distribution is empty, has a zero or non-finite
    /// total weight, omits the lightest declared isotope, contains an undeclared
    /// mass number or is not strictly increasing in nominal mass. The generator
    /// is then unchanged.
    pub fn set_isotope_override(
        &mut self,
        symbol: &str,
        distribution: IsotopeDistribution,
    ) -> Result<()> {
        let mut position = 0;
        let atom = parse_atom(symbol.as_bytes(), &mut position)?;
        if position != symbol.len() {
            return Err(invalid("expected one element or isotope symbol"));
        }
        let natural = atom_distribution(atom);
        let dense = dense_from_distribution(&distribution, None)?;
        let natural_dense = dense_from_distribution(&natural, None)?;
        if distribution.is_empty()
            || distribution.probability_sum() <= 0.0
            || !distribution.probability_sum().is_finite()
        {
            return Err(invalid(
                "isotope override needs finite, positive total probability",
            ));
        }
        if dense.origin != natural_dense.origin {
            return Err(invalid(
                "override must include the lightest declared isotope (zero weight is allowed)",
            ));
        }
        let allowed: BTreeSet<_> = natural
            .peaks
            .iter()
            .map(|p| p.mass.round() as u64)
            .collect();
        if distribution
            .peaks
            .iter()
            .any(|p| !allowed.contains(&(p.mass.round() as u64)))
        {
            return Err(invalid(
                "override contains an undeclared isotope mass number",
            ));
        }
        self.overrides.insert(atom, distribution);
        Ok(())
    }

    /// Remove every isotope override (source `clearIsotopeOverrides`).
    pub fn clear_isotope_overrides(&mut self) {
        self.overrides.clear();
    }

    /// Generate a normalized pattern. Rejects negative atom counts/charge.
    ///
    /// Positive charge follows the pinned C++ behavior: convolve extra natural H
    /// atoms, then use the proton-mass shift in the lightest-mass anchor. Returned
    /// coordinates are masses, never divided by charge. For an explicit neutral
    /// adduct formula, add its atom counts and set charge to zero before calling.
    ///
    /// Source `run` iterates through all elements, convolves each element's
    /// distribution according to its atom count, replaces nominal bins by
    /// corrected masses and renormalizes. The source marks the implicit charge
    /// handling deprecated: OpenMS 4 is announced to ignore the charge, but the
    /// pinned source still adds `charge` hydrogen atoms, and so does this port.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative charge (source
    /// `Exception::Precondition`) or atom count, when a vector or the work
    /// budget exceeds its limit, a probability overflows, or the retained
    /// probabilities sum to zero.
    pub fn run(&self, formula: &EmpiricalFormula) -> Result<IsotopeDistribution> {
        self.run_with_work(formula, &mut CoarseIsotopeWork::default())
    }

    /// Share probability-product accounting across theoretical envelopes.
    pub(crate) fn run_with_work(
        &self,
        formula: &EmpiricalFormula,
        work: &mut CoarseIsotopeWork,
    ) -> Result<IsotopeDistribution> {
        self.check_run_preconditions(formula)?;
        if self.precision == ProbabilityPrecision::SourceSingle {
            return self.run_source_single(formula, work);
        }
        let mut pattern = DensePattern::identity();
        for (&atom, &count) in &formula.atoms {
            let distribution = self.distribution_for(atom);
            let atomic = dense_from_distribution(&distribution, self.max_peaks)?;
            let repeated = power_dense(&atomic, count as usize, self.max_peaks, work)?;
            pattern = convolve_dense(&pattern, &repeated, self.max_peaks, work)?;
        }
        if formula.charge > 0 {
            let hydrogen = Atom {
                symbol: "H",
                isotope: None,
            };
            let distribution = self.distribution_for(hydrogen);
            let atomic = dense_from_distribution(&distribution, self.max_peaks)?;
            let repeated = power_dense(&atomic, formula.charge as usize, self.max_peaks, work)?;
            pattern = convolve_dense(&pattern, &repeated, self.max_peaks, work)?;
        }
        let mut result = self.correct_masses(&pattern.probabilities, lightest_mass(formula))?;
        result.renormalize()?;
        Ok(result)
    }

    /// `CoarseIsotopePatternGenerator.cpp:76-132` in binary32.
    ///
    /// The source convolves `H^charge` even for charge zero; that convolution
    /// with the identity copies every binary32 value exactly and is skipped.
    fn run_source_single(
        &self,
        formula: &EmpiricalFormula,
        work: &mut CoarseIsotopeWork,
    ) -> Result<IsotopeDistribution> {
        let mut result = self.run_source_single_unnormalized(formula, work)?;
        result.renormalize_with(ProbabilityPrecision::SourceSingle)?;
        Ok(result)
    }

    /// [`Self::estimate_from_peptide_weight`] under
    /// [`ProbabilityPrecision::SourceSingle`] for a caller that follows the
    /// source past a zero probability sum.
    ///
    /// Where every retained binary32 bin underflows to zero, source
    /// `renormalize` divides each by the zero sum
    /// (`IsotopeDistribution.cpp:176-192`), so every weight is NaN; this
    /// returns [`SourceSingleEstimate::AllUnderflowed`] there instead of the
    /// error [`Self::estimate_from_peptide_weight`] returns. Everything else is
    /// that function's result. Only `FeatureFinderAlgorithmPicked` step 2.5
    /// calls it; the other callers keep the error.
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_from_peptide_weight`], except for the all-underflow
    /// case, and [`Error::InvalidValue`] when the generator does not use
    /// [`ProbabilityPrecision::SourceSingle`].
    pub(crate) fn estimate_from_peptide_weight_source(
        &self,
        average_mass: f64,
    ) -> Result<SourceSingleEstimate> {
        if self.precision != ProbabilityPrecision::SourceSingle {
            return Err(invalid(
                "the source estimate requires ProbabilityPrecision::SourceSingle",
            ));
        }
        let formula = AveragineComposition::PEPTIDE
            .estimate_average_mass(average_mass)?
            .formula;
        let mut work = CoarseIsotopeWork::default();
        self.check_run_preconditions(&formula)?;
        let mut result = self.run_source_single_unnormalized(&formula, &mut work)?;
        let mut all_zero = true;
        for peak in result.peaks() {
            if narrow(peak.probability)? != 0.0 {
                all_zero = false;
            }
        }
        if all_zero && !result.is_empty() {
            return Ok(SourceSingleEstimate::AllUnderflowed { len: result.len() });
        }
        result.renormalize_with(ProbabilityPrecision::SourceSingle)?;
        Ok(SourceSingleEstimate::Normalized(result))
    }

    /// The checks [`Self::run`] applies before any convolution.
    fn check_run_preconditions(&self, formula: &EmpiricalFormula) -> Result<()> {
        if formula.charge < 0 {
            return Err(invalid(
                "coarse isotope generation does not support negative charge",
            ));
        }
        if formula.atoms.values().any(|&count| count < 0) {
            return Err(invalid(
                "isotope generation requires nonnegative atom counts",
            ));
        }
        Ok(())
    }

    /// [`Self::run_source_single`] before its final `renormalize`.
    fn run_source_single_unnormalized(
        &self,
        formula: &EmpiricalFormula,
        work: &mut CoarseIsotopeWork,
    ) -> Result<IsotopeDistribution> {
        let atoms = source_atom_order(formula);
        let mut pattern = SinglePattern::identity();
        for &(atom, count) in &atoms {
            let isotopes = single_from_distribution(&self.distribution_for(atom), self.max_peaks)?;
            let repeated = single_power(&isotopes, count as usize, self.max_peaks, work)?;
            pattern = single_convolve(&pattern, &repeated, self.max_peaks, work)?;
        }
        if formula.charge > 0 {
            let hydrogen = Atom {
                symbol: "H",
                isotope: None,
            };
            let isotopes =
                single_from_distribution(&self.distribution_for(hydrogen), self.max_peaks)?;
            let repeated = single_power(&isotopes, formula.charge as usize, self.max_peaks, work)?;
            pattern = single_convolve(&pattern, &repeated, self.max_peaks, work)?;
        }
        let lightest = source_lightest_mass(&atoms, formula.charge);
        let probabilities: Vec<f64> = pattern.intensities.iter().copied().map(f64::from).collect();
        self.correct_masses(&probabilities, lightest)
    }

    fn distribution_for(&self, atom: Atom) -> IsotopeDistribution {
        self.overrides
            .get(&atom)
            .cloned()
            .unwrap_or_else(|| atom_distribution(atom))
    }

    /// Lightest-isotope mass anchor in the order of the selected precision.
    fn lightest_mass_for(&self, formula: &EmpiricalFormula) -> f64 {
        match self.precision {
            ProbabilityPrecision::Double => lightest_mass(formula),
            ProbabilityPrecision::SourceSingle => {
                source_lightest_mass(&source_atom_order(formula), formula.charge)
            }
        }
    }

    /// Nominal convolution: round input masses, fill gaps with zero weights, and
    /// truncate to max_peaks. Inputs must be strictly increasing after rounding.
    /// Output weights are unnormalized and output masses are integer bin sums.
    ///
    /// Source public `convolve` (`CoarseIsotopePatternGenerator.cpp:316-357`)
    /// accumulates from the highest indices down so small products come first.
    /// Under [`ProbabilityPrecision::SourceSingle`] each accumulation step is the
    /// source binary32 `p + l * r`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when input nominal masses are not
    /// strictly increasing or exceed the exact integer range, a vector or the
    /// work budget exceeds its limit, a weight exceeds the binary32 range in
    /// source precision, or a probability overflows.
    pub fn convolve(
        &self,
        left: &IsotopeDistribution,
        right: &IsotopeDistribution,
    ) -> Result<IsotopeDistribution> {
        if self.precision == ProbabilityPrecision::SourceSingle {
            let left = single_from_distribution(left, self.max_peaks)?;
            let right = single_from_distribution(right, self.max_peaks)?;
            return single_convolve(
                &left,
                &right,
                self.max_peaks,
                &mut CoarseIsotopeWork::default(),
            )?
            .into_distribution();
        }
        let left = dense_from_distribution(left, self.max_peaks)?;
        let right = dense_from_distribution(right, self.max_peaks)?;
        let result = convolve_dense(
            &left,
            &right,
            self.max_peaks,
            &mut CoarseIsotopeWork::default(),
        )?;
        result.into_distribution()
    }

    /// Repeated nominal convolution by exponentiation by squaring.
    /// A nonempty input raised to zero gives the identity; empty stays empty.
    ///
    /// Source protected `convolvePow_` (`CoarseIsotopePatternGenerator.cpp:359-419`)
    /// returns its input unchanged for exponent one; this returns the gap-filled
    /// nominal form for every exponent.
    ///
    /// # Errors
    ///
    /// As [`Self::convolve`].
    pub fn convolve_power(
        &self,
        input: &IsotopeDistribution,
        exponent: usize,
    ) -> Result<IsotopeDistribution> {
        if self.precision == ProbabilityPrecision::SourceSingle {
            let input = single_from_distribution(input, self.max_peaks)?;
            return single_power(
                &input,
                exponent,
                self.max_peaks,
                &mut CoarseIsotopeWork::default(),
            )?
            .into_distribution();
        }
        let input = dense_from_distribution(input, self.max_peaks)?;
        power_dense(
            &input,
            exponent,
            self.max_peaks,
            &mut CoarseIsotopeWork::default(),
        )?
        .into_distribution()
    }

    fn correct_masses(&self, probabilities: &[f64], lightest: f64) -> Result<IsotopeDistribution> {
        validate_nonnegative(lightest, "lightest isotope mass")?;
        let peaks = probabilities
            .iter()
            .enumerate()
            .map(|(i, &probability)| {
                let mass = lightest + i as f64 * C13C12_MASSDIFF_U;
                IsotopePeak {
                    mass: if self.mass_mode == CoarseMassMode::Nominal {
                        mass.round()
                    } else {
                        mass
                    },
                    probability,
                }
            })
            .collect();
        IsotopeDistribution::from_peaks(peaks)
    }

    /// Estimate a peptide distribution from its average weight (source
    /// `estimateFromPeptideWeight`).
    ///
    /// Uses the averagine model of Senko et al., "Determination of Monoisotopic
    /// Masses and Ion Populations for Large Biomolecules from Resolved Isotopic
    /// Distributions", as tabulated in [`AveragineComposition::PEPTIDE`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative or non-finite weight, an
    /// atom count outside `i32`, or any error of [`Self::run`].
    pub fn estimate_from_peptide_weight(&self, average_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_weight_and_comp(average_mass, AveragineComposition::PEPTIDE)
    }

    /// Estimate a peptide distribution from its monoisotopic weight (source
    /// `estimateFromPeptideMonoWeight`).
    ///
    /// Same Senko averagine model, but the formula is fitted to the monoisotopic
    /// mass, so no monoisotopic mass determination is performed.
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_from_peptide_weight`].
    pub fn estimate_from_peptide_mono_weight(&self, mono_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_mono_weight_and_comp(mono_mass, AveragineComposition::PEPTIDE)
    }

    /// Estimate a nucleotide distribution from its average weight (source
    /// `estimateFromRNAWeight`).
    ///
    /// Uses the averagine model of Zubarev and Demirev, "Isotope depletion of
    /// large biomolecules: Implications for molecular mass measurements", as
    /// tabulated in [`AveragineComposition::RNA`].
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_from_peptide_weight`].
    pub fn estimate_from_rna_weight(&self, average_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_weight_and_comp(average_mass, AveragineComposition::RNA)
    }

    /// Estimate a nucleotide distribution from its monoisotopic weight (source
    /// `estimateFromRNAMonoWeight`), with the Zubarev and Demirev RNA averagine.
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_from_peptide_weight`].
    pub fn estimate_from_rna_mono_weight(&self, mono_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_mono_weight_and_comp(mono_mass, AveragineComposition::RNA)
    }

    /// Estimate a nucleotide distribution from its average weight (source
    /// `estimateFromDNAWeight`), with the Zubarev and Demirev DNA averagine in
    /// [`AveragineComposition::DNA`].
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_from_peptide_weight`].
    pub fn estimate_from_dna_weight(&self, average_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_weight_and_comp(average_mass, AveragineComposition::DNA)
    }

    /// Estimate a distribution from an average weight and an average composition
    /// (source `estimateFromWeightAndComp`).
    ///
    /// The formula comes from [`AveragineComposition::estimate_average_mass`]:
    /// rounded heavy-element counts, then hydrogen fitted to the remaining mass.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative or non-finite weight or
    /// relative count, a composition without positive mass, an atom count
    /// outside `i32`, or any error of [`Self::run`].
    pub fn estimate_from_weight_and_comp(
        &self,
        average_mass: f64,
        composition: AveragineComposition,
    ) -> Result<IsotopeDistribution> {
        self.run(&composition.estimate_average_mass(average_mass)?.formula)
    }

    /// Estimate a distribution from a monoisotopic weight and an average
    /// composition (source `estimateFromMonoWeightAndComp`).
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_from_weight_and_comp`].
    pub fn estimate_from_mono_weight_and_comp(
        &self,
        mono_mass: f64,
        composition: AveragineComposition,
    ) -> Result<IsotopeDistribution> {
        self.run(&composition.estimate_mono_mass(mono_mass)?.formula)
    }

    /// Estimate a peptide distribution from its average weight and an exact
    /// sulfur count (source `estimateFromPeptideWeightAndS`).
    ///
    /// The remaining mass uses the sulfur-free Senko averagine. The source
    /// preconditions are `sulfur_count <= average_mass / average_weight(S)` and
    /// `average_mass >= 0`; both are checked here.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the weight is negative or not
    /// finite, the sulfur mass exceeds the weight, or any error of
    /// [`Self::run`] occurs.
    pub fn estimate_from_peptide_weight_and_sulfur(
        &self,
        average_mass: f64,
        sulfur_count: u32,
    ) -> Result<IsotopeDistribution> {
        self.run(
            &AveragineComposition::PEPTIDE
                .estimate_average_mass_with_sulfur(average_mass, sulfur_count)?
                .formula,
        )
    }

    /// Joint fragment/isolation weights using Rockwood's conditional-isotope rule.
    ///
    /// Inputs must enumerate successive isotope indices starting at M0 (including
    /// zero bins); masses are ignored. `precursor_isotopes` are zero-based indices
    /// and duplicates are ignored. Normalize the result to get conditional
    /// probabilities. Truncation uses the output bound throughout, avoiding the
    /// upstream out-of-bounds loop when max_peaks is smaller than the fragment.
    ///
    /// Source `calcFragmentIsotopeDist` follows Rockwood, Kushnir and Nelson,
    /// "Dissociation of Individual Isotopic Peaks: Predicting Isotopic
    /// Distributions of Product Ions in MSn". Its precondition that both inputs
    /// are gapless is the M0-indexed layout required above.
    /// `fragment_lightest_mass` anchors the output masses; the source names it
    /// the fragment's monoisotopic mass. Under
    /// [`ProbabilityPrecision::SourceSingle`] the complementary weights are
    /// accumulated in binary32 in ascending precursor index and the sum is
    /// multiplied by the binary32 fragment weight.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when no precursor isotope is selected, an
    /// index or the output exceeds its limit, the lightest mass is negative or
    /// not finite, or a weight exceeds the binary32 range in source precision.
    pub fn calc_fragment_isotope_dist(
        &self,
        fragment: &IsotopeDistribution,
        complementary: &IsotopeDistribution,
        precursor_isotopes: &[usize],
        fragment_lightest_mass: f64,
    ) -> Result<IsotopeDistribution> {
        validate_nonnegative(fragment_lightest_mass, "fragment lightest isotope mass")?;
        let selected = selected_isotopes(precursor_isotopes)?;
        if fragment.is_empty() || complementary.is_empty() {
            return Ok(IsotopeDistribution::empty());
        }
        let length = bounded_length(fragment.len(), self.max_peaks)?;
        let products = length
            .checked_mul(selected.len())
            .ok_or_else(|| invalid("fragment work size overflow"))?;
        CoarseIsotopeWork::default().consume(products)?;
        let mut probabilities = Vec::with_capacity(length);
        if self.precision == ProbabilityPrecision::SourceSingle {
            // CoarseIsotopePatternGenerator.cpp:509-520, Peak1D float slots.
            for index in 0..length {
                let weight = selected
                    .iter()
                    .filter_map(|&precursor| precursor.checked_sub(index))
                    .filter_map(|index| complementary.peaks.get(index))
                    .try_fold(0.0_f32, |sum, peak| {
                        narrow(peak.probability).map(|value| sum + value)
                    })?;
                let joint = weight * narrow(fragment.peaks[index].probability)?;
                if !joint.is_finite() {
                    return Err(invalid("isotope probability overflow in fragment weights"));
                }
                probabilities.push(f64::from(joint));
            }
            return self.correct_masses(&probabilities, fragment_lightest_mass);
        }
        for index in 0..length {
            let weight: f64 = selected
                .iter()
                .filter_map(|&precursor| precursor.checked_sub(index))
                .filter_map(|index| complementary.peaks.get(index))
                .map(|p| p.probability)
                .sum();
            probabilities.push(fragment.peaks[index].probability * weight);
        }
        self.correct_masses(&probabilities, fragment_lightest_mass)
    }

    /// Estimate fragment and complementary formulas independently from average
    /// weights, then return unnormalized isolation-conditioned fragment weights.
    ///
    /// Covers source `estimateForFragmentFromPeptideWeight`,
    /// `estimateForFragmentFromRNAWeight`, `estimateForFragmentFromDNAWeight`
    /// and `estimateForFragmentFromWeightAndComp` through `composition`. As in
    /// the source, the inner solver keeps this generator's isotope overrides and
    /// precision and reports `max(precursor_isotopes) + 1` bins; the output mass
    /// anchor is the fragment formula's lightest-isotope mass. The source
    /// preconditions (both weights positive, fragment not heavier than the
    /// precursor, at least one isotope) are checked.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when a weight is not positive and finite,
    /// the fragment is heavier than the precursor, no isotope is selected, or
    /// any error of [`Self::run`] or [`Self::calc_fragment_isotope_dist`]
    /// occurs.
    pub fn estimate_fragment_from_weights(
        &self,
        precursor_mass: f64,
        fragment_mass: f64,
        precursor_isotopes: &[usize],
        composition: AveragineComposition,
    ) -> Result<IsotopeDistribution> {
        validate_fragment_weights(precursor_mass, fragment_mass)?;
        let selected = selected_isotopes(precursor_isotopes)?;
        let mut solver = self.clone();
        solver.max_peaks = Some(selected.last().copied().expect("nonempty") + 1);
        let formula = composition.estimate_average_mass(fragment_mass)?.formula;
        let fragment = solver.run(&formula)?;
        let complementary =
            solver.estimate_from_weight_and_comp(precursor_mass - fragment_mass, composition)?;
        self.calc_fragment_isotope_dist(
            &fragment,
            &complementary,
            precursor_isotopes,
            self.lightest_mass_for(&formula),
        )
    }

    /// Poisson approximation from Bellew et al., with lambda = mass / 1800.
    /// Retains the source recurrence and normalization order for finite weights;
    /// falls back to log space when the direct weights or their sum overflow.
    ///
    /// Source `approximateIntensities` (Bellew et al.,
    /// <https://dx.doi.org/10.1093/bioinformatics/btl276>) is documented as
    /// roughly 100 times faster than [`Self::estimate_from_peptide_weight`] and
    /// only an approximation; `num_peaks` is independent of the generator's
    /// maximum isotope count. This computes in `f64`, as the source does; it is
    /// not affected by [`ProbabilityPrecision`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative or non-finite mass, zero
    /// peaks, or more than [`MAX_ISOTOPE_PEAKS`] peaks.
    pub fn approximate_intensities(mass: f64, num_peaks: usize) -> Result<Vec<f64>> {
        validate_nonnegative(mass, "peptide mass")?;
        if num_peaks == 0 {
            return Err(invalid("Poisson approximation requires at least one peak"));
        }
        check_size(num_peaks)?;
        if mass == 0.0 {
            let mut result = vec![0.0; num_peaks];
            result[0] = 1.0;
            return Ok(result);
        }
        let factor = mass / 1800.0;
        let mut probabilities = vec![1.0; num_peaks];
        let mut current = 1.0;
        let mut sum = 1.0;
        for (index, probability) in probabilities.iter_mut().enumerate().skip(1) {
            // Division precedes multiplication in the source; its last bits
            // matter when the model is used at a KL acceptance threshold.
            current *= factor / index as f64;
            *probability = current;
            sum += current;
            if !sum.is_finite() {
                break;
            }
        }
        if sum.is_finite() {
            for probability in &mut probabilities {
                *probability /= sum;
            }
            return Ok(probabilities);
        }
        let log_factor = mass.ln() - 1800.0_f64.ln();
        probabilities.fill(0.0);
        for index in 1..num_peaks {
            probabilities[index] = probabilities[index - 1] + log_factor - (index as f64).ln();
        }
        let maximum = probabilities
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        for probability in &mut probabilities {
            *probability = (*probability - maximum).exp();
        }
        let sum: f64 = probabilities.iter().sum();
        for probability in &mut probabilities {
            *probability /= sum;
        }
        Ok(probabilities)
    }

    /// Upstream Poisson mass layout: first mass is exactly `mass`, with neutron
    /// spacing / charge; mass itself is not divided by charge or protonated.
    ///
    /// Source `approximateFromPeptideWeight` is documented as about 50 times
    /// faster than [`Self::estimate_from_peptide_weight`]. For monoisotopic mass
    /// 1000 the source lists first intensities 0.573753, 0.318752, 0.0885422
    /// against 0.571133, 0.306181, 0.0958111 from the averagine estimate; KL
    /// divergences over 20 peaks range from 4.97e-5 at mass 20 to 0.0144 at mass
    /// 2500, below the 0.05 (two peaks) and 0.6 (six or more peaks) isotope
    /// pattern thresholds of Teo et al. The source keeps a `float` running
    /// product; this uses the `f64` recurrence of
    /// [`Self::approximate_intensities`] and is not affected by
    /// [`ProbabilityPrecision`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for zero charge or any error of
    /// [`Self::approximate_intensities`].
    pub fn approximate_from_peptide_weight(
        mass: f64,
        num_peaks: usize,
        charge: u32,
    ) -> Result<IsotopeDistribution> {
        if charge == 0 {
            return Err(invalid("Poisson mass spacing requires positive charge"));
        }
        let probabilities = Self::approximate_intensities(mass, num_peaks)?;
        IsotopeDistribution::from_peaks(
            probabilities
                .into_iter()
                .enumerate()
                .map(|(i, probability)| IsotopePeak {
                    mass: mass + i as f64 * NEUTRON_MASS_U / f64::from(charge),
                    probability,
                })
                .collect(),
        )
    }
}

/// Relative elemental counts used for average-composition mass estimates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AveragineComposition {
    /// Relative carbon count.
    pub carbon: f64,
    /// Relative hydrogen count.
    pub hydrogen: f64,
    /// Relative nitrogen count.
    pub nitrogen: f64,
    /// Relative oxygen count.
    pub oxygen: f64,
    /// Relative sulfur count.
    pub sulfur: f64,
    /// Relative phosphorus count.
    pub phosphorus: f64,
}

/// Rounded heavy-element formula plus hydrogen mass adjustment.
#[derive(Clone, Debug, PartialEq)]
pub struct FormulaEstimate {
    /// The estimated formula, charge zero.
    pub formula: EmpiricalFormula,
    /// False when rounded heavy elements exceed the requested mass enough to
    /// need negative hydrogen; the returned formula then contains no hydrogen.
    pub hydrogen_adjustment_succeeded: bool,
}

impl AveragineComposition {
    /// Senko peptide averagine, as tabulated by OpenMS.
    pub const PEPTIDE: Self = Self {
        carbon: 4.9384,
        hydrogen: 7.7583,
        nitrogen: 1.3577,
        oxygen: 1.4773,
        sulfur: 0.0417,
        phosphorus: 0.0,
    };
    /// Zubarev and Demirev nucleotide averagine for RNA, as tabulated by OpenMS.
    pub const RNA: Self = Self {
        carbon: 9.75,
        hydrogen: 12.25,
        nitrogen: 3.75,
        oxygen: 7.0,
        sulfur: 0.0,
        phosphorus: 1.0,
    };
    /// Zubarev and Demirev nucleotide averagine for DNA (one oxygen fewer than
    /// RNA), as tabulated by OpenMS.
    pub const DNA: Self = Self {
        carbon: 9.75,
        hydrogen: 12.25,
        nitrogen: 3.75,
        oxygen: 6.0,
        sulfur: 0.0,
        phosphorus: 1.0,
    };

    /// Estimate a formula whose average mass approximates `mass` (source
    /// `EmpiricalFormula::estimateFromWeightAndComp`).
    ///
    /// Heavy-element counts are the relative counts scaled to the mass and
    /// rounded half away from zero; hydrogen is then fitted to the remaining
    /// mass. The source returns `false` when that needs negative hydrogen;
    /// this reports it in [`FormulaEstimate::hydrogen_adjustment_succeeded`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a negative or non-finite mass or
    /// relative count, a composition without positive unit mass, or a count
    /// outside `i32`.
    pub fn estimate_average_mass(self, mass: f64) -> Result<FormulaEstimate> {
        self.estimate(mass, true)
    }
    /// Estimate a formula whose monoisotopic mass approximates `mass` (source
    /// `EmpiricalFormula::estimateFromMonoWeightAndComp`).
    ///
    /// # Errors
    ///
    /// As [`Self::estimate_average_mass`].
    pub fn estimate_mono_mass(self, mass: f64) -> Result<FormulaEstimate> {
        self.estimate(mass, false)
    }

    /// Fix sulfur count, estimating the remaining average mass with sulfur-free
    /// composition. This safely inserts S even when its estimated count is zero.
    ///
    /// Source `EmpiricalFormula::estimateFromWeightAndCompAndS`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the mass is negative or not finite,
    /// the sulfur count overflows `i32`, the sulfur mass exceeds `mass`, or any
    /// error of [`Self::estimate_average_mass`] occurs.
    pub fn estimate_average_mass_with_sulfur(
        mut self,
        mass: f64,
        sulfur_count: u32,
    ) -> Result<FormulaEstimate> {
        validate_nonnegative(mass, "target mass")?;
        let sulfur_count =
            i32::try_from(sulfur_count).map_err(|_| invalid("sulfur count overflows i32"))?;
        let remaining =
            mass - f64::from(sulfur_count) * element("S").expect("table sulfur").average_mass();
        if remaining < 0.0 {
            return Err(invalid("fixed sulfur mass exceeds target mass"));
        }
        self.sulfur = 0.0;
        let mut result = self.estimate_average_mass(remaining)?;
        if sulfur_count > 0 {
            result.formula.atoms.insert(
                Atom {
                    symbol: "S",
                    isotope: None,
                },
                sulfur_count,
            );
        }
        Ok(result)
    }

    fn estimate(self, mass: f64, average: bool) -> Result<FormulaEstimate> {
        validate_nonnegative(mass, "target mass")?;
        let components = [
            ("C", self.carbon),
            ("H", self.hydrogen),
            ("N", self.nitrogen),
            ("O", self.oxygen),
            ("S", self.sulfur),
            ("P", self.phosphorus),
        ];
        let mut unit_mass = 0.0;
        for (symbol, relative) in components {
            validate_nonnegative(relative, "averagine relative count")?;
            let entry = element(symbol).expect("averagine table element");
            unit_mass += relative
                * if average {
                    entry.average_mass()
                } else {
                    entry.mono_mass()
                };
        }
        if !unit_mass.is_finite() || unit_mass <= 0.0 {
            return Err(invalid("averagine composition needs finite, positive mass"));
        }
        let scale = mass / unit_mass;
        let mut formula = EmpiricalFormula::default();
        for (symbol, relative) in components {
            if symbol == "H" {
                continue;
            }
            let count = rounded_count(relative * scale)?;
            if count > 0 {
                formula.atoms.insert(
                    Atom {
                        symbol,
                        isotope: None,
                    },
                    count,
                );
            }
        }
        let hydrogen_mass = element("H").expect("table hydrogen");
        let residual = mass
            - if average {
                formula.average_mass()
            } else {
                formula.mono_mass()
            };
        let hydrogen = (residual
            / if average {
                hydrogen_mass.average_mass()
            } else {
                hydrogen_mass.mono_mass()
            })
        .round();
        let succeeded = hydrogen >= 0.0;
        if succeeded {
            let count = rounded_count(hydrogen)?;
            if count > 0 {
                formula.atoms.insert(
                    Atom {
                        symbol: "H",
                        isotope: None,
                    },
                    count,
                );
            }
        }
        Ok(FormulaEstimate {
            formula,
            hydrogen_adjustment_succeeded: succeeded,
        })
    }
}

fn rounded_count(value: f64) -> Result<i32> {
    let rounded = value.round();
    if !rounded.is_finite() || rounded < 0.0 || rounded > f64::from(i32::MAX) {
        return Err(invalid("estimated atom count is outside i32 range"));
    }
    Ok(rounded as i32)
}

fn atom_distribution(atom: Atom) -> IsotopeDistribution {
    let entry = element(atom.symbol).expect("validated formula element");
    let peaks = entry
        .isotopes()
        .iter()
        .filter(|i| atom.isotope.is_none_or(|number| i.mass_number == number))
        .map(|isotope| IsotopePeak {
            mass: isotope.mass,
            probability: if atom.isotope.is_some() {
                1.0
            } else {
                isotope.abundance
            },
        })
        .collect();
    IsotopeDistribution { peaks }
}

fn lightest_mass(formula: &EmpiricalFormula) -> f64 {
    formula.atoms.iter().fold(
        f64::from(formula.charge) * PROTON_MASS_U,
        |mass, (&atom, &count)| mass + f64::from(count) * lightest_isotope_mass(atom),
    )
}

/// The first tabulated isotope mass of an element, or the labelled isotope's mass.
fn lightest_isotope_mass(atom: Atom) -> f64 {
    let entry = element(atom.symbol).expect("validated formula element");
    match atom.isotope {
        Some(number) => {
            entry
                .isotopes()
                .iter()
                .find(|i| i.mass_number == number)
                .expect("validated isotope")
                .mass
        }
        None => entry.isotopes()[0].mass,
    }
}

/// Formula atoms in the source-precision convolution order: ascending atomic
/// number, each natural element before its labelled isotopes (ascending mass
/// number).
///
/// The source's `std::map<const Element*, _>` iterates in `ElementDB` address
/// order, which varies between runs of the same binary. Ascending atomic number
/// is the majority order measured for natural elements; the labelled-isotope
/// placement is a native choice, since no placement was a majority.
fn source_atom_order(formula: &EmpiricalFormula) -> Vec<(Atom, i32)> {
    let mut atoms: Vec<(Atom, i32)> = formula
        .atoms
        .iter()
        .map(|(&atom, &count)| (atom, count))
        .collect();
    atoms.sort_by_key(|&(atom, _)| {
        let entry = element(atom.symbol).expect("validated formula element");
        (entry.atomic_number, atom.isotope.unwrap_or(0))
    });
    atoms
}

/// `EmpiricalFormula::getLightestIsotopeWeight` (`EmpiricalFormula.cpp:57-67`)
/// summed in source element order.
fn source_lightest_mass(atoms: &[(Atom, i32)], charge: i32) -> f64 {
    atoms
        .iter()
        .fold(PROTON_MASS_U * f64::from(charge), |mass, &(atom, count)| {
            mass + lightest_isotope_mass(atom) * f64::from(count)
        })
}

#[derive(Clone, Debug)]
struct DensePattern {
    origin: f64,
    probabilities: Vec<f64>,
}

impl DensePattern {
    fn identity() -> Self {
        Self {
            origin: 0.0,
            probabilities: vec![1.0],
        }
    }
    fn empty() -> Self {
        Self {
            origin: 0.0,
            probabilities: Vec::new(),
        }
    }
    fn into_distribution(self) -> Result<IsotopeDistribution> {
        IsotopeDistribution::from_peaks(
            self.probabilities
                .into_iter()
                .enumerate()
                .map(|(i, probability)| IsotopePeak {
                    mass: self.origin + i as f64,
                    probability,
                })
                .collect(),
        )
    }
}

fn dense_from_distribution(
    distribution: &IsotopeDistribution,
    max_peaks: Option<usize>,
) -> Result<DensePattern> {
    if distribution.is_empty() {
        return Ok(DensePattern::empty());
    }
    let origin = distribution.peaks[0].mass.round();
    let mut previous = None;
    for peak in &distribution.peaks {
        let nominal = peak.mass.round();
        if nominal > MAX_EXACT_INTEGER {
            return Err(invalid("nominal isotope mass exceeds exact integer range"));
        }
        if previous.is_some_and(|last| nominal <= last) {
            return Err(invalid(
                "convolution requires strictly increasing, distinct nominal masses",
            ));
        }
        previous = Some(nominal);
    }
    let span = previous.expect("nonempty distribution") - origin;
    if span > MAX_ISOTOPE_PEAKS as f64 && max_peaks.is_none() {
        return Err(invalid("isotope gap span exceeds allocation limit"));
    }
    let length = bounded_length((span as usize).saturating_add(1), max_peaks)?;
    let mut probabilities = vec![0.0; length];
    for peak in &distribution.peaks {
        let index = (peak.mass.round() - origin) as usize;
        if index < length {
            probabilities[index] = peak.probability;
        }
    }
    Ok(DensePattern {
        origin,
        probabilities,
    })
}

fn convolve_dense(
    left: &DensePattern,
    right: &DensePattern,
    max_peaks: Option<usize>,
    work: &mut CoarseIsotopeWork,
) -> Result<DensePattern> {
    if left.probabilities.is_empty() || right.probabilities.is_empty() {
        return Ok(DensePattern::empty());
    }
    let full_length = left
        .probabilities
        .len()
        .checked_add(right.probabilities.len())
        .and_then(|n| n.checked_sub(1))
        .ok_or_else(|| invalid("isotope convolution length overflow"))?;
    let length = bounded_length(full_length, max_peaks)?;
    let origin = left.origin + right.origin;
    if origin + length as f64 > MAX_EXACT_INTEGER {
        return Err(invalid(
            "nominal convolution mass exceeds exact integer range",
        ));
    }
    let left_length = left.probabilities.len().min(length);
    let mut products = 0_usize;
    for i in 0..left_length {
        products = products
            .checked_add(right.probabilities.len().min(length - i))
            .ok_or_else(|| invalid("convolution work size overflow"))?;
    }
    work.consume(products)?;
    let mut probabilities = vec![0.0; length];
    for i in (0..left_length).rev() {
        for j in (0..right.probabilities.len().min(length - i)).rev() {
            probabilities[i + j] += left.probabilities[i] * right.probabilities[j];
            if !probabilities[i + j].is_finite() {
                return Err(invalid("isotope probability overflow in convolution"));
            }
        }
    }
    Ok(DensePattern {
        origin,
        probabilities,
    })
}

fn power_dense(
    input: &DensePattern,
    mut exponent: usize,
    max_peaks: Option<usize>,
    work: &mut CoarseIsotopeWork,
) -> Result<DensePattern> {
    if input.probabilities.is_empty() {
        return Ok(DensePattern::empty());
    }
    let mut result = DensePattern::identity();
    let mut power = input.clone();
    while exponent > 0 {
        if exponent & 1 != 0 {
            result = convolve_dense(&result, &power, max_peaks, work)?;
        }
        exponent >>= 1;
        if exponent > 0 {
            power = convolve_dense(&power, &power, max_peaks, work)?;
        }
    }
    Ok(result)
}

/// Binary32 coarse bins laid out like the source's gap-filled `Peak1D`
/// containers: a nominal origin mass and one intensity per consecutive bin.
#[derive(Clone, Debug)]
struct SinglePattern {
    origin: f64,
    intensities: Vec<f32>,
}

impl SinglePattern {
    fn identity() -> Self {
        Self {
            origin: 0.0,
            intensities: vec![1.0],
        }
    }
    fn empty() -> Self {
        Self {
            origin: 0.0,
            intensities: Vec::new(),
        }
    }
    fn into_distribution(self) -> Result<IsotopeDistribution> {
        IsotopeDistribution::from_peaks(
            self.intensities
                .into_iter()
                .enumerate()
                .map(|(i, intensity)| IsotopePeak {
                    mass: self.origin + i as f64,
                    probability: f64::from(intensity),
                })
                .collect(),
        )
    }
}

/// Narrow a validated weight as a `float` intensity assignment does.
fn narrow(value: f64) -> Result<f32> {
    let single = value as f32;
    if !single.is_finite() {
        return Err(invalid(
            "isotope probability exceeds the binary32 range of the source Peak1D intensity",
        ));
    }
    Ok(single)
}

/// Source `fillGaps_` (`CoarseIsotopePatternGenerator.cpp:525-547`) followed by
/// the `float` intensity narrowing, truncated to `max_peaks` bins.
///
/// The source does not truncate inputs, but a bin at or beyond `max_peaks` never
/// contributes to a retained bin (see [`single_convolve`]), so the retained
/// binary32 values are unchanged.
fn single_from_distribution(
    distribution: &IsotopeDistribution,
    max_peaks: Option<usize>,
) -> Result<SinglePattern> {
    let dense = dense_from_distribution(distribution, max_peaks)?;
    let intensities = dense
        .probabilities
        .iter()
        .map(|&probability| narrow(probability))
        .collect::<Result<Vec<f32>>>()?;
    Ok(SinglePattern {
        origin: dense.origin,
        intensities,
    })
}

/// Source `convolve` and `convolveSquare_` (`CoarseIsotopePatternGenerator.cpp:
/// 316-357, 421-448`) in binary32.
///
/// Bins are accumulated for `i` descending and, within each `i`, `j` descending,
/// each step `p + l[i] * r[j]` as a separate `f32` product and sum. The result
/// keeps at most `max_peaks` bins. The source squares keep `max_peaks + 1` bins;
/// retained bin `k < max_peaks` receives the products `(i, k - i)` for `i` from
/// `k` down to 0 in both layouts, and the extra bin only feeds bins at or beyond
/// `max_peaks`, so every retained value is identical.
fn single_convolve(
    left: &SinglePattern,
    right: &SinglePattern,
    max_peaks: Option<usize>,
    work: &mut CoarseIsotopeWork,
) -> Result<SinglePattern> {
    if left.intensities.is_empty() || right.intensities.is_empty() {
        return Ok(SinglePattern::empty());
    }
    let full_length = left
        .intensities
        .len()
        .checked_add(right.intensities.len())
        .and_then(|n| n.checked_sub(1))
        .ok_or_else(|| invalid("isotope convolution length overflow"))?;
    let length = bounded_length(full_length, max_peaks)?;
    let origin = left.origin + right.origin;
    if origin + length as f64 > MAX_EXACT_INTEGER {
        return Err(invalid(
            "nominal convolution mass exceeds exact integer range",
        ));
    }
    let left_length = left.intensities.len().min(length);
    let mut products = 0_usize;
    for i in 0..left_length {
        products = products
            .checked_add(right.intensities.len().min(length - i))
            .ok_or_else(|| invalid("convolution work size overflow"))?;
    }
    work.consume(products)?;
    let mut intensities = vec![0.0_f32; length];
    for i in (0..left_length).rev() {
        for j in (0..right.intensities.len().min(length - i)).rev() {
            let product = left.intensities[i] * right.intensities[j];
            let sum = intensities[i + j] + product;
            if !sum.is_finite() {
                return Err(invalid("isotope probability overflow in convolution"));
            }
            intensities[i + j] = sum;
        }
    }
    Ok(SinglePattern {
        origin,
        intensities,
    })
}

/// Source `convolvePow_` (`CoarseIsotopePatternGenerator.cpp:359-419`) in binary32.
///
/// Exponent one returns the input. Otherwise the result starts as the input
/// (odd exponent) or the identity, and is convolved with each successive square
/// whose bit is set, from the lowest bit up. Squares the source computes after
/// the highest set bit are never used and are not computed.
fn single_power(
    input: &SinglePattern,
    exponent: usize,
    max_peaks: Option<usize>,
    work: &mut CoarseIsotopeWork,
) -> Result<SinglePattern> {
    if input.intensities.is_empty() {
        return Ok(SinglePattern::empty());
    }
    if exponent == 1 {
        return Ok(input.clone());
    }
    let mut result = if exponent & 1 == 1 {
        input.clone()
    } else {
        SinglePattern::identity()
    };
    let mut square = input.clone();
    let mut remaining = exponent >> 1;
    while remaining > 0 {
        square = single_convolve(&square, &square, max_peaks, work)?;
        if remaining & 1 == 1 {
            result = single_convolve(&result, &square, max_peaks, work)?;
        }
        remaining >>= 1;
    }
    Ok(result)
}

/// Probability-product budget shared by the convolutions of one calculation.
#[derive(Default)]
pub(crate) struct CoarseIsotopeWork {
    used: usize,
}
impl CoarseIsotopeWork {
    /// Charge `products` against [`MAX_CONVOLUTION_PRODUCTS`].
    pub(crate) fn consume(&mut self, products: usize) -> Result<()> {
        self.used = self
            .used
            .checked_add(products)
            .ok_or_else(|| invalid("isotope work count overflow"))?;
        if self.used > MAX_CONVOLUTION_PRODUCTS {
            return Err(invalid(
                "isotope convolution exceeds work limit; request fewer peaks",
            ));
        }
        Ok(())
    }
}

fn bounded_length(length: usize, max_peaks: Option<usize>) -> Result<usize> {
    let length = max_peaks.map_or(length, |maximum| length.min(maximum));
    check_size(length)?;
    Ok(length)
}

fn check_size(length: usize) -> Result<()> {
    if length > MAX_ISOTOPE_PEAKS {
        return Err(invalid("isotope vector exceeds allocation limit"));
    }
    Ok(())
}

fn selected_isotopes(indices: &[usize]) -> Result<BTreeSet<usize>> {
    if indices.is_empty() {
        return Err(invalid("at least one precursor isotope must be selected"));
    }
    if indices.len() > MAX_ISOTOPE_PEAKS || indices.iter().any(|&i| i >= MAX_ISOTOPE_PEAKS) {
        return Err(invalid(
            "selected precursor isotope exceeds work/allocation limit",
        ));
    }
    Ok(indices.iter().copied().collect())
}

fn validate_fragment_weights(precursor: f64, fragment: f64) -> Result<()> {
    validate_nonnegative(precursor, "precursor mass")?;
    validate_nonnegative(fragment, "fragment mass")?;
    if precursor <= 0.0 || fragment <= 0.0 || fragment > precursor {
        return Err(invalid(
            "fragment estimation requires 0 < fragment mass <= precursor mass",
        ));
    }
    Ok(())
}

fn validate_nonnegative(value: f64, name: &str) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(invalid(format!("{name} must be finite and nonnegative")));
    }
    Ok(())
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn products(precision: ProbabilityPrecision, max_peaks: Option<usize>, text: &str) -> usize {
        let generator = CoarseIsotopePatternGenerator::new(max_peaks, CoarseMassMode::Approximate)
            .unwrap()
            .with_precision(precision);
        let formula: EmpiricalFormula = text.parse().unwrap();
        let mut work = CoarseIsotopeWork::default();
        generator.run_with_work(&formula, &mut work).unwrap();
        work.used
    }

    /// Source precision charges the products of its own convolution sequence.
    /// It copies instead of convolving with the identity for odd exponents, and
    /// it convolves in atomic-number rather than symbol order, which changes the
    /// intermediate lengths. Its total is therefore lower or higher than the
    /// native total depending on the formula and the bin limit.
    #[test]
    fn source_precision_work_follows_its_own_convolution_sequence() {
        use ProbabilityPrecision::{Double, SourceSingle};
        // Exponent one: natively C^1 costs 2 (identity * C) plus 2 (pattern * C^1);
        // source precision copies C for the power and pays only the second 2.
        assert_eq!(products(Double, None, "C1"), 4);
        assert_eq!(products(SourceSingle, None, "C1"), 2);
        // Unbounded peptide: H^95 (191 bins) is convolved before C^44 (45 bins),
        // which costs more than the skipped identity convolutions save.
        assert_eq!(products(Double, None, "C44H95N12O13S1"), 36_245);
        assert_eq!(products(SourceSingle, None, "C44H95N12O13S1"), 36_380);
        // Twenty bins cap both lengths, so only the skipped convolutions remain.
        assert_eq!(products(Double, Some(20), "C44H95N12O13S1"), 3_081);
        assert_eq!(products(SourceSingle, Some(20), "C44H95N12O13S1"), 3_070);
    }
}
