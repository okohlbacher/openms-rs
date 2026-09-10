// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Internal b/a fragments retain the source's start, length and loss conventions.

use super::*;

impl TheoreticalSpectrumGenerator {
    pub(super) fn internal_fragment_peaks(
        &self,
        rows: &mut Vec<PeakRow>,
        peptide: &AASequence,
        min_charge: u8,
        max_charge: u8,
        remaining_loss_entries: &mut usize,
    ) -> Result<()> {
        let water = EmpiricalFormula::parse("H2O")?;
        let carbon_monoxide = EmpiricalFormula::parse("CO")?;
        for charge in min_charge..=max_charge {
            for series in [TheoreticalIonSeries::B, TheoreticalIonSeries::A] {
                let is_a = series == TheoreticalIonSeries::A;
                let intensity = self.intensities.for_series(series);
                let offset = if is_a {
                    -carbon_monoxide.mono_mass()
                } else {
                    0.0
                };
                // Source l+3<n deliberately omits the last possible dipeptide.
                for start in 1..peptide.len().saturating_sub(3) {
                    let end = (start + 10).min(peptide.len() - 1);
                    let mut mass = PROTON_MASS_U * f64::from(charge);
                    let mut ladder = Vec::with_capacity(end - start - 1);
                    for index in start..end {
                        mass += peptide.internal_mass(index)?.ok_or_else(|| {
                            Error::Unsupported("internal fragment has an unresolved mass".into())
                        })?;
                        if index == start {
                            continue;
                        }
                        let fragment = peptide.subsequence(start..index + 1)?;
                        let mut formula = match fragment.formula() {
                            Ok(formula) => Some(formula.checked_sub(&water)?),
                            Err(Error::Unsupported(_)) => None,
                            Err(error) => return Err(error),
                        };
                        if is_a {
                            formula = formula
                                .map(|f| f.checked_sub(&carbon_monoxide))
                                .transpose()?;
                        }
                        if formula.as_ref().is_some_and(|f| !physical_formula(f)) {
                            return Err(invalid("internal ion produces a negative atom count"));
                        }
                        let name =
                            format!("{}{}", fragment.as_str(), if is_a { "-CO" } else { "" });
                        let losses = if self.add_losses {
                            // The first residue's loss declarations are not visited
                            // by the source loop; terminal losses are never added.
                            fragment_losses(&fragment, series, false, 1, remaining_loss_entries)?
                        } else {
                            Vec::new()
                        };
                        push_peak(
                            rows,
                            (mass + offset) / f64::from(charge),
                            intensity,
                            &name,
                            charge,
                        )?;
                        ladder.push(IonTemplate {
                            name,
                            formula,
                            // Unlike terminal templates, this includes protons:
                            // preserve source addition order before loss subtraction.
                            mono_mass: mass + offset,
                            losses,
                        });
                    }
                    // Each start emits all intact lengths before their loss peaks.
                    for ion in ladder {
                        for loss in ion.losses {
                            if ion
                                .formula
                                .as_ref()
                                .map(|f| f.checked_sub(&loss))
                                .transpose()?
                                .is_some_and(|f| !physical_formula(&f))
                            {
                                continue;
                            }
                            push_peak(
                                rows,
                                (ion.mono_mass - loss.mono_mass()) / f64::from(charge),
                                intensity * self.relative_loss_intensity,
                                &format!(
                                    "{}-{}{}",
                                    ion.name,
                                    loss,
                                    "+".repeat(usize::from(charge))
                                ),
                                charge,
                            )?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
