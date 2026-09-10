// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Coarse isotope distributions and averagine estimates.
//!
//! Ported from OpenMS4-core revision `7c029e8`, principally
//! `ISOTOPEDISTRIBUTION/IsotopeDistribution.cpp`,
//! `ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.cpp`, and
//! `CHEMISTRY/EmpiricalFormula.cpp`. Probabilities use f64 rather than C++'s f32.
//! Coarse mass correction is a carbon-13 spacing approximation, not isotope fine
//! structure. See `docs/ISOTOPE_SUPPORT.md` for conventions and limits.

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
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsotopePeak {
    pub mass: f64,
    pub probability: f64,
}

/// Validated isotope masses and weights, preserving insertion order.
///
/// The default is the convolution identity `(0, 1)`. Use [`Self::empty`] for an
/// empty container. Mutating methods preserve finite, nonnegative weights and
/// masses; read-only peak slices cannot bypass validation.
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
    pub fn empty() -> Self {
        Self { peaks: Vec::new() }
    }

    pub fn from_peaks(peaks: Vec<IsotopePeak>) -> Result<Self> {
        check_size(peaks.len())?;
        for peak in &peaks {
            validate_nonnegative(peak.mass, "isotope mass")?;
            validate_nonnegative(peak.probability, "isotope probability")?;
        }
        Ok(Self { peaks })
    }

    pub fn peaks(&self) -> &[IsotopePeak] {
        &self.peaks
    }
    pub fn len(&self) -> usize {
        self.peaks.len()
    }
    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty()
    }
    pub fn clear(&mut self) {
        self.peaks.clear();
    }

    pub fn insert(&mut self, peak: IsotopePeak) -> Result<()> {
        check_size(self.len().saturating_add(1))?;
        validate_nonnegative(peak.mass, "isotope mass")?;
        validate_nonnegative(peak.probability, "isotope probability")?;
        self.peaks.push(peak);
        Ok(())
    }

    /// Resize, padding with `(0, 0)` as in OpenMS. Does not preserve sortedness.
    pub fn resize(&mut self, length: usize) -> Result<()> {
        check_size(length)?;
        self.peaks.resize(length, IsotopePeak::default());
        Ok(())
    }

    pub fn min_mass(&self) -> Option<f64> {
        self.peaks.iter().map(|p| p.mass).min_by(f64::total_cmp)
    }

    pub fn max_mass(&self) -> Option<f64> {
        self.peaks.iter().map(|p| p.mass).max_by(f64::total_cmp)
    }

    /// Most abundant peak; the first wins ties. Empty distributions return None.
    pub fn most_abundant(&self) -> Option<IsotopePeak> {
        self.peaks.iter().copied().reduce(|best, peak| {
            if peak.probability > best.probability {
                peak
            } else {
                best
            }
        })
    }

    pub fn sort_by_mass(&mut self) {
        self.peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
    }

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

    /// Probability-weighted mass, without requiring prior normalization.
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
    pub fn trim_intensities(&mut self, cutoff: f64) -> Result<()> {
        validate_nonnegative(cutoff, "probability cutoff")?;
        self.peaks.retain(|p| p.probability >= cutoff);
        Ok(())
    }

    /// Remove the low-probability prefix in current order, without normalization.
    /// Unlike upstream's all-below-cutoff bug, this can remove every peak.
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

    /// Remove the low-probability suffix in current order, without normalization.
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

/// Mass labels for the same nominal-isotope probabilities.
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
/// `run` normalizes retained probabilities to one. `convolve`, `convolve_power`,
/// and `calc_fragment_isotope_dist` do not normalize. None means all peaks,
/// subject to documented allocation/work limits; Some(0) is invalid.
#[derive(Clone, Debug, Default)]
pub struct CoarseIsotopePatternGenerator {
    max_peaks: Option<usize>,
    mass_mode: CoarseMassMode,
    overrides: BTreeMap<Atom, IsotopeDistribution>,
}

impl CoarseIsotopePatternGenerator {
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
            overrides: BTreeMap::new(),
        })
    }

    pub fn max_peaks(&self) -> Option<usize> {
        self.max_peaks
    }
    pub fn mass_mode(&self) -> CoarseMassMode {
        self.mass_mode
    }

    /// Override isotope weights locally, without mutating the shared element table.
    /// Masses may be exact or nominal but must round to declared isotope numbers.
    /// Include the lightest declared isotope, even at zero probability, so the
    /// mass origin remains defined. Weights need not sum to one.
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

    pub fn clear_isotope_overrides(&mut self) {
        self.overrides.clear();
    }

    /// Generate a normalized pattern. Rejects negative atom counts/charge.
    ///
    /// Positive charge follows the pinned C++ behavior: convolve extra natural H
    /// atoms, then use the proton-mass shift in the lightest-mass anchor. Returned
    /// coordinates are masses, never divided by charge. For an explicit neutral
    /// adduct formula, add its atom counts and set charge to zero before calling.
    pub fn run(&self, formula: &EmpiricalFormula) -> Result<IsotopeDistribution> {
        self.run_with_work(formula, &mut CoarseIsotopeWork::default())
    }

    /// Share probability-product accounting across theoretical envelopes.
    pub(crate) fn run_with_work(
        &self,
        formula: &EmpiricalFormula,
        work: &mut CoarseIsotopeWork,
    ) -> Result<IsotopeDistribution> {
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

    fn distribution_for(&self, atom: Atom) -> IsotopeDistribution {
        self.overrides
            .get(&atom)
            .cloned()
            .unwrap_or_else(|| atom_distribution(atom))
    }

    /// Nominal convolution: round input masses, fill gaps with zero weights, and
    /// truncate to max_peaks. Inputs must be strictly increasing after rounding.
    /// Output weights are unnormalized and output masses are integer bin sums.
    pub fn convolve(
        &self,
        left: &IsotopeDistribution,
        right: &IsotopeDistribution,
    ) -> Result<IsotopeDistribution> {
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
    pub fn convolve_power(
        &self,
        input: &IsotopeDistribution,
        exponent: usize,
    ) -> Result<IsotopeDistribution> {
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

    pub fn estimate_from_peptide_weight(&self, average_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_weight_and_comp(average_mass, AveragineComposition::PEPTIDE)
    }

    pub fn estimate_from_peptide_mono_weight(&self, mono_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_mono_weight_and_comp(mono_mass, AveragineComposition::PEPTIDE)
    }

    pub fn estimate_from_rna_weight(&self, average_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_weight_and_comp(average_mass, AveragineComposition::RNA)
    }

    pub fn estimate_from_rna_mono_weight(&self, mono_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_mono_weight_and_comp(mono_mass, AveragineComposition::RNA)
    }

    pub fn estimate_from_dna_weight(&self, average_mass: f64) -> Result<IsotopeDistribution> {
        self.estimate_from_weight_and_comp(average_mass, AveragineComposition::DNA)
    }

    pub fn estimate_from_weight_and_comp(
        &self,
        average_mass: f64,
        composition: AveragineComposition,
    ) -> Result<IsotopeDistribution> {
        self.run(&composition.estimate_average_mass(average_mass)?.formula)
    }

    pub fn estimate_from_mono_weight_and_comp(
        &self,
        mono_mass: f64,
        composition: AveragineComposition,
    ) -> Result<IsotopeDistribution> {
        self.run(&composition.estimate_mono_mass(mono_mass)?.formula)
    }

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
            lightest_mass(&formula),
        )
    }

    /// Poisson approximation from Bellew et al., with lambda = mass / 1800.
    /// Retains the source recurrence and normalization order for finite weights;
    /// falls back to log space when the direct weights or their sum overflow.
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
    pub carbon: f64,
    pub hydrogen: f64,
    pub nitrogen: f64,
    pub oxygen: f64,
    pub sulfur: f64,
    pub phosphorus: f64,
}

/// Rounded heavy-element formula plus hydrogen mass adjustment.
#[derive(Clone, Debug, PartialEq)]
pub struct FormulaEstimate {
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
    pub const RNA: Self = Self {
        carbon: 9.75,
        hydrogen: 12.25,
        nitrogen: 3.75,
        oxygen: 7.0,
        sulfur: 0.0,
        phosphorus: 1.0,
    };
    pub const DNA: Self = Self {
        carbon: 9.75,
        hydrogen: 12.25,
        nitrogen: 3.75,
        oxygen: 6.0,
        sulfur: 0.0,
        phosphorus: 1.0,
    };

    pub fn estimate_average_mass(self, mass: f64) -> Result<FormulaEstimate> {
        self.estimate(mass, true)
    }
    pub fn estimate_mono_mass(self, mass: f64) -> Result<FormulaEstimate> {
        self.estimate(mass, false)
    }

    /// Fix sulfur count, estimating the remaining average mass with sulfur-free
    /// composition. This safely inserts S even when its estimated count is zero.
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
        |mass, (&atom, &count)| {
            let entry = element(atom.symbol).expect("validated formula element");
            let lightest = match atom.isotope {
                Some(number) => {
                    entry
                        .isotopes()
                        .iter()
                        .find(|i| i.mass_number == number)
                        .expect("validated isotope")
                        .mass
                }
                None => entry.isotopes()[0].mass,
            };
            mass + f64::from(count) * lightest
        },
    )
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

#[derive(Default)]
pub(crate) struct CoarseIsotopeWork {
    used: usize,
}
impl CoarseIsotopeWork {
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
