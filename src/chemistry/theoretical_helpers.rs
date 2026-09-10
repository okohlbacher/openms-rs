// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Compact m/z-only helper from TheoreticalSpectrumGenerator.cpp:1165–1261
//! at OpenMS4-core revision 7c029e8. Kept separate from intensity/annotation work.

use super::{
    AASequence, EmpiricalFormula, MAX_RESIDUE_WORK, MAX_THEORETICAL_PEAKS,
    MAX_THEORETICAL_RESIDUES, PROTON_MASS_U, TheoreticalIonSeries, TheoreticalSpectrumGenerator,
    invalid,
};
use crate::{Error, Result};

impl TheoreticalSpectrumGenerator {
    /// Append compact m/z values for all enabled ordinary a/b/c/x/y/z series,
    /// charges `charge` down through one and lengths one through the FULL peptide.
    /// Each ladder uses its corresponding terminal modification only. The
    /// existing and new values are always sorted together, including when the
    /// peptide is empty or `charge` is zero (both generate no new values).
    ///
    /// This is the source `getPrefixAndSuffixIonsMZ` helper. Ion-series flags are
    /// its only spectrum settings: radical z variants, first-prefix selection,
    /// losses, intensities, precursors, isotopes, annotations and sort options are
    /// ignored. Repeated series act as a single enabled flag. Monoisotopic mass
    /// tags need no invented composition; unresolved required masses are errors.
    ///
    /// Limits are 100,000 combined values, 4,096 residues when generating ions
    /// and 10,000,000 estimated scan/ladder/sort work units. Existing/new values
    /// must be finite and nonnegative, including after conversion to f32. Every
    /// fallible step precedes appending, so errors leave output unchanged.
    pub fn append_mass_spectrum(
        &self,
        output: &mut Vec<f32>,
        peptide: &AASequence,
        charge: u8,
    ) -> Result<()> {
        if output.len() > MAX_THEORETICAL_PEAKS || self.ion_series.len() > MAX_RESIDUE_WORK {
            return Err(invalid(
                "mass spectrum exceeds output or settings work limit",
            ));
        }
        // Source bool flags, in its fixed charge-descending helper order.
        let ordered = [
            TheoreticalIonSeries::B,
            TheoreticalIonSeries::Y,
            TheoreticalIonSeries::A,
            TheoreticalIonSeries::X,
            TheoreticalIonSeries::C,
            TheoreticalIonSeries::Z,
        ];
        let mut enabled = [false; 6];
        for series in &self.ion_series {
            match series {
                TheoreticalIonSeries::B => enabled[0] = true,
                TheoreticalIonSeries::Y => enabled[1] = true,
                TheoreticalIonSeries::A => enabled[2] = true,
                TheoreticalIonSeries::X => enabled[3] = true,
                TheoreticalIonSeries::C => enabled[4] = true,
                TheoreticalIonSeries::Z => enabled[5] = true,
                TheoreticalIonSeries::ZPlusOne | TheoreticalIonSeries::ZPlusTwo => {}
            }
        }
        let series_count = enabled.iter().filter(|&&value| value).count();
        let active = charge != 0 && series_count != 0 && !peptide.is_empty();
        if active && peptide.len() > MAX_THEORETICAL_RESIDUES {
            return Err(invalid("peptide exceeds mass spectrum residue limit"));
        }
        let additions = if active {
            peptide
                .len()
                .checked_mul(series_count)
                .and_then(|value| value.checked_mul(usize::from(charge)))
                .ok_or_else(|| invalid("mass spectrum peak count overflow"))?
        } else {
            0
        };
        let combined = output
            .len()
            .checked_add(additions)
            .filter(|&count| count <= MAX_THEORETICAL_PEAKS)
            .ok_or_else(|| invalid("mass spectrum exceeds combined peak limit"))?;
        let sort_depth = usize::BITS - combined.leading_zeros();
        let work = combined
            .checked_mul(sort_depth as usize)
            .and_then(|value| value.checked_add(self.ion_series.len()))
            .and_then(|value| value.checked_add(additions))
            .and_then(|value| value.checked_add(if active { peptide.len() } else { 0 }))
            .ok_or_else(|| invalid("mass spectrum work overflow"))?;
        if work > MAX_RESIDUE_WORK {
            return Err(invalid("mass spectrum exceeds work limit"));
        }
        if output
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "existing mass spectrum values must be finite and nonnegative",
            ));
        }
        let mut added = Vec::with_capacity(additions);
        if active {
            // Compute once using the same validated residue mass conventions as
            // sequence/fragment chemistry. Never subtract a large prefix from a
            // full peptide to recover a small suffix.
            let masses: Vec<f64> = (0..peptide.len())
                .map(|index| {
                    peptide.internal_mass(index)?.ok_or_else(|| {
                        Error::Unsupported(
                            "mass spectrum requires a known mass at every residue".into(),
                        )
                    })
                })
                .collect::<Result<_>>()?;
            let water = EmpiricalFormula::parse("H2O")?;
            for z in (1..=charge).rev() {
                for (series, enabled) in ordered.into_iter().zip(enabled) {
                    if !enabled {
                        continue;
                    }
                    let terminal = if series.is_prefix() {
                        peptide.n_terminal_modification()
                    } else {
                        peptide.c_terminal_modification()
                    };
                    // Preserve the source's f64 operation order until the final
                    // f32 store: protons, terminal delta, ion conversion, residues.
                    let mut mass = PROTON_MASS_U * f64::from(z);
                    if let Some(modification) = terminal {
                        mass += modification.diff_mono_mass()?;
                    }
                    mass += series.formula_delta().checked_add(&water)?.mono_mass();
                    for ordinal in 0..masses.len() {
                        let index = if series.is_prefix() {
                            ordinal
                        } else {
                            masses.len() - ordinal - 1
                        };
                        mass += masses[index];
                        let mz = mass / f64::from(z);
                        let stored = mz as f32;
                        if !mz.is_finite() || mz < 0.0 || !stored.is_finite() {
                            return Err(invalid(
                                "generated mass spectrum value is not finite and nonnegative f32",
                            ));
                        }
                        added.push(stored);
                    }
                }
            }
        }
        output.extend(added);
        output.sort_by(f32::total_cmp);
        Ok(())
    }
}
