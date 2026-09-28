// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The remainder of `CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h`
//! on [`CoarseIsotopePatternGenerator`]: the two setters, the override map
//! accessor and the sulfur-aware fragment estimate.
//!
//! # API mapping
//!
//! | Source (`CoarseIsotopePatternGenerator.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `setMaxIsotope(const Size&)` | [`CoarseIsotopePatternGenerator::set_max_peaks`] |
//! | `setRoundMasses(bool)` | [`CoarseIsotopePatternGenerator::set_mass_mode`] |
//! | `getIsotopeOverrides()` | [`CoarseIsotopePatternGenerator::isotope_overrides`] |
//! | `estimateForFragmentFromPeptideWeightAndS` | [`CoarseIsotopePatternGenerator::estimate_fragment_from_peptide_weight_and_sulfur`] |

use super::{
    AveragineComposition, CoarseIsotopePatternGenerator, CoarseMassMode, IsotopeDistribution,
    check_size, invalid, selected_isotopes, validate_fragment_weights,
};
use crate::Result;

impl CoarseIsotopePatternGenerator {
    /// Replaces the maximum number of reported isotopes: source
    /// `setMaxIsotope(const Size&)`, where the source's `0` ("all isotopes") is
    /// `None`.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] for `Some(0)` and for a count
    /// above [`super::MAX_ISOTOPE_PEAKS`], as [`Self::new`] does; the generator
    /// is then unchanged.
    pub fn set_max_peaks(&mut self, max_peaks: Option<usize>) -> Result<()> {
        if let Some(length) = max_peaks {
            if length == 0 {
                return Err(invalid("maximum isotope count must be positive"));
            }
            check_size(length)?;
        }
        self.max_peaks = max_peaks;
        Ok(())
    }

    /// Replaces how output masses are labelled: source `setRoundMasses(bool)`,
    /// `true` being [`CoarseMassMode::Nominal`].
    pub fn set_mass_mode(&mut self, mass_mode: CoarseMassMode) {
        self.mass_mode = mass_mode;
    }

    /// The isotope overrides set with [`Self::set_isotope_override`], as
    /// (symbol, distribution) pairs: source `getIsotopeOverrides()`.
    ///
    /// The source map is keyed by `ElementDB` element pointer and iterates in
    /// address order; this yields the element or isotope symbol (`C`, `(13)C`,
    /// and `(2)H` for an override set as `D`) in the fixed order of the
    /// formula entries (symbol, natural before labelled, mass number).
    pub fn isotope_overrides(
        &self,
    ) -> impl ExactSizeIterator<Item = (String, &IsotopeDistribution)> {
        self.overrides.iter().map(|(atom, distribution)| {
            let symbol = match atom.isotope {
                Some(number) => format!("({number}){}", atom.symbol),
                None => atom.symbol.to_owned(),
            };
            (symbol, distribution)
        })
    }

    /// Isolation-conditioned fragment weights from peptide average weights and
    /// sulfur counts: source `estimateForFragmentFromPeptideWeightAndS`.
    ///
    /// Following `CoarseIsotopePatternGenerator.cpp:197-218`, the inner solver
    /// reports `max(precursor_isotopes) + 1` bins in approximate mass mode and
    /// keeps this generator's isotope overrides (and here its precision); the
    /// fragment formula is the peptide averagine estimate with
    /// `sulfur_fragment` sulfur atoms, the complementary fragment is
    /// [`Self::estimate_from_peptide_weight_and_sulfur`] of the remaining weight
    /// and sulfur, and [`Self::calc_fragment_isotope_dist`] combines them,
    /// anchored at the fragment formula's lightest-isotope mass. The result is
    /// not normalized.
    ///
    /// # Errors
    ///
    /// Returns [`crate::Error::InvalidValue`] when a weight is not positive and
    /// finite, the fragment is heavier than the precursor, no isotope is
    /// selected, `sulfur_fragment` exceeds `sulfur_precursor` (the source
    /// subtracts them as unsigned integers and wraps), or any error of the
    /// estimates, [`Self::run`] or [`Self::calc_fragment_isotope_dist`] occurs.
    pub fn estimate_fragment_from_peptide_weight_and_sulfur(
        &self,
        precursor_mass: f64,
        sulfur_precursor: u32,
        fragment_mass: f64,
        sulfur_fragment: u32,
        precursor_isotopes: &[usize],
    ) -> Result<IsotopeDistribution> {
        validate_fragment_weights(precursor_mass, fragment_mass)?;
        let selected = selected_isotopes(precursor_isotopes)?;
        let sulfur_complementary = sulfur_precursor
            .checked_sub(sulfur_fragment)
            .ok_or_else(|| invalid("fragment sulfur count exceeds the precursor's"))?;
        let deepest = selected
            .last()
            .copied()
            .ok_or_else(|| invalid("at least one precursor isotope must be selected"))?;
        let mut solver = self.clone();
        solver.max_peaks = Some(deepest + 1);
        solver.mass_mode = CoarseMassMode::Approximate;
        let formula = AveragineComposition::PEPTIDE
            .estimate_average_mass_with_sulfur(fragment_mass, sulfur_fragment)?
            .formula;
        let fragment = solver.run(&formula)?;
        let complementary = solver.estimate_from_peptide_weight_and_sulfur(
            precursor_mass - fragment_mass,
            sulfur_complementary,
        )?;
        self.calc_fragment_isotope_dist(
            &fragment,
            &complementary,
            precursor_isotopes,
            self.lightest_mass_for(&formula),
        )
    }
}
