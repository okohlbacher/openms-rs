// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Conservative preflight for the exact monoisotopic ProForma wrapper settings.
use super::*;
use crate::chemistry::SequenceModification;
use std::mem::size_of;

struct Budget<'a> {
    work: &'a mut usize,
    bytes: &'a mut usize,
}
impl Budget<'_> {
    fn charge(&mut self, work: usize, bytes: usize) -> Result<()> {
        *self.work = self.work.checked_sub(work).ok_or_else(limit)?;
        *self.bytes = self.bytes.checked_sub(bytes).ok_or_else(limit)?;
        Ok(())
    }
}
fn limit() -> Error {
    invalid("ProForma theoretical generation resource limit exceeded")
}
fn tree(entries: usize) -> usize {
    if entries == 0 {
        0
    } else {
        512usize.saturating_add(
            entries
                .saturating_mul(3)
                .max(11)
                .saturating_mul(size_of::<(crate::chemistry::Atom, i32)>().saturating_add(32)),
        )
    }
}
impl TheoreticalSpectrumGenerator {
    /// This deliberately narrow adapter is for source ProForma's selected flags.
    /// Standalone generation and its category limits retain their existing behavior.
    pub(crate) fn generate_for_proforma(
        &self,
        peptide: &AASequence,
        min: u8,
        max: u8,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<MSSpectrum> {
        if self.isotope_model != TheoreticalIsotopeModel::None
            || self.add_internal_fragments
            || self.add_terminal_losses
            || self.add_all_precursor_charges
        {
            return Err(invalid(
                "unsupported options for private ProForma generation budget",
            ));
        }
        let mut budget = Budget { work, bytes };
        let n = peptide.len();
        if n > MAX_THEORETICAL_RESIDUES || self.ion_series.len() > 6 {
            return Err(limit());
        }
        budget.charge(n.saturating_add(128), 4096)?;
        // Sequence slices share Known handles; measure owned spelling/formula
        // state before the existing helper walks annotation slots.
        let payload = peptide.generation_payload_bytes()?;
        budget.charge(payload, payload)?;
        let mut atoms = 6usize;
        let mut losses = 0usize;
        let mut loss_bytes = 0usize;
        let mut loss_atoms = 0usize;
        let mut largest_loss = 0usize;
        for index in 0..n {
            let modification = peptide.residue_modification(index)?;
            if let Some(SequenceModification::Known(record)) = modification {
                atoms = atoms
                    .saturating_add(record.diff_formula().stored_atom_types())
                    .saturating_add(
                        record
                            .absolute_formula()
                            .map_or(0, EmpiricalFormula::stored_atom_types),
                    );
                if self.add_losses {
                    let count = record.neutral_losses().len();
                    budget.charge(count, 0)?;
                    losses = losses.saturating_add(count);
                    for loss in record.neutral_losses() {
                        let count = loss.formula().stored_atom_types();
                        budget.charge(count, 0)?;
                        loss_atoms = loss_atoms.saturating_add(count);
                        largest_loss = largest_loss.max(count);
                        loss_bytes = loss_bytes
                            .saturating_add(tree(count).saturating_mul(3))
                            .saturating_add(count.saturating_add(1).saturating_mul(128))
                            .saturating_add(2048);
                    }
                }
            } else if modification.is_none() && self.add_losses {
                let count = match peptide.as_str().as_bytes()[index] {
                    b'R' => 3,
                    b'D' | b'E' | b'S' | b'T' | b'K' | b'N' | b'Q' => 1,
                    _ => 0,
                };
                losses = losses.saturating_add(count);
                loss_atoms = loss_atoms.saturating_add(count * 3);
                largest_loss = largest_loss.max(3);
                loss_bytes = loss_bytes.saturating_add(count.saturating_mul(tree(3) * 3 + 4096));
            }
        }
        for modification in [
            peptide.n_terminal_modification(),
            peptide.c_terminal_modification(),
        ]
        .into_iter()
        .flatten()
        {
            if let Some(record) = modification.known() {
                atoms = atoms
                    .saturating_add(record.diff_formula().stored_atom_types())
                    .saturating_add(
                        record
                            .absolute_formula()
                            .map_or(0, EmpiricalFormula::stored_atom_types),
                    );
            }
        }
        let map_bytes = tree(atoms);
        let templates = n.saturating_sub(1).saturating_mul(self.ion_series.len());
        let charges = usize::from(max.saturating_sub(min)) + 1;
        // Every fragment reconstructs chemistry from retained residues. Use the
        // complete peptide's map/payload bound even for shorter slices. Cover
        // per-residue internal formula maps, accumulated maps, mass queries,
        // template maps/vectors, and duplicate loss-key formatting/comparisons.
        let levels = (usize::BITS - losses.max(1).leading_zeros()) as usize;
        let name = 32usize
            .saturating_add(usize::from(max))
            .saturating_add(largest_loss.saturating_add(1).saturating_mul(64));
        if name > 4 * 1024 * 1024 {
            return Err(limit());
        }
        let rebuild_work = n
            .saturating_mul(atoms.saturating_mul(16).saturating_add(256))
            .saturating_add(payload);
        let rebuild_bytes = n
            .saturating_mul(map_bytes.saturating_mul(8).saturating_add(1024))
            .saturating_add(payload);
        let loss_work = losses
            .saturating_mul(levels.saturating_add(1))
            .saturating_mul(name.saturating_mul(4))
            .saturating_add(loss_atoms.saturating_mul(32));
        budget.charge(
            templates
                .saturating_add(4)
                .saturating_mul(rebuild_work.saturating_add(loss_work)),
            templates.saturating_add(4).saturating_mul(
                rebuild_bytes
                    .saturating_add(loss_bytes)
                    .saturating_add(4096),
            ),
        )?;
        // Loss formula subtraction and formatting run per charge even when a
        // negative atom count suppresses the resulting peak.
        let loss_visits = templates.saturating_mul(losses).saturating_mul(charges);
        budget.charge(
            loss_visits.saturating_mul(
                atoms
                    .saturating_add(largest_loss)
                    .saturating_mul(32)
                    .saturating_add(name),
            ),
            loss_visits.saturating_mul(
                tree(atoms.saturating_add(largest_loss)).saturating_add(name.saturating_mul(2)),
            ),
        )?;
        let rows = templates
            .saturating_mul(losses.saturating_add(1))
            .saturating_mul(charges)
            .saturating_add(usize::from(self.add_precursor_peaks) * 3)
            .saturating_add(usize::from(self.add_abundant_immonium_ions) * 7)
            .min(MAX_THEORETICAL_PEAKS);
        let levels = (usize::BITS - rows.max(1).leading_zeros()) as usize;
        budget.charge(
            rows.saturating_mul(levels.saturating_add(name).saturating_add(32)),
            rows.saturating_mul(
                size_of::<PeakRow>()
                    .saturating_mul(3)
                    .saturating_add(size_of::<Peak1D>() + size_of::<String>() + size_of::<i32>())
                    .saturating_add(name.saturating_mul(4)),
            )
            .saturating_add(32768),
        )?;
        self.generate(peptide, min, max, None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn failure_after_preflight_keeps_consumed_counters_and_public_behavior() {
        let generator = TheoreticalSpectrumGenerator::default();
        let sequence = AASequence::parse("PEPTIDE").unwrap();
        let expected = generator.generate(&sequence, 1, 2, None).unwrap();
        let mut work = 50_000_000;
        let mut bytes = 256 * 1024 * 1024;
        assert_eq!(
            generator
                .generate_for_proforma(&sequence, 1, 2, &mut work, &mut bytes)
                .unwrap(),
            expected
        );
        let used = [50_000_000 - work, 256 * 1024 * 1024 - bytes];
        for index in 0..2 {
            let mut work = if index == 0 { used[0] - 1 } else { 50_000_000 };
            let mut bytes = if index == 1 {
                used[1] - 1
            } else {
                256 * 1024 * 1024
            };
            let before = (work, bytes);
            assert!(
                generator
                    .generate_for_proforma(&sequence, 1, 2, &mut work, &mut bytes)
                    .is_err()
            );
            assert!(work < before.0 && bytes < before.1);
        }
        let mut work = 50_000_000;
        let mut bytes = 256 * 1024 * 1024;
        assert!(
            generator
                .generate_for_proforma(&sequence, 0, 1, &mut work, &mut bytes)
                .is_err()
        );
        assert!(work < 50_000_000 && bytes < 256 * 1024 * 1024);
        assert_eq!(generator.generate(&sequence, 1, 2, None).unwrap(), expected);
    }
    #[test]
    fn duplicate_custom_loss_declarations_are_charged_before_generation() {
        use crate::chemistry::{
            ModificationRecord, ModificationsDB, NeutralLoss, ResidueModification,
        };
        let formula = EmpiricalFormula::parse("H2O").unwrap();
        let loss =
            NeutralLoss::new(formula.clone(), formula.mono_mass(), formula.average_mass()).unwrap();
        let registry = ModificationsDB::from_records(vec![
            ResidueModification::from_record(ModificationRecord {
                name: "Bulk".into(),
                origin: Some('K'),
                neutral_losses: vec![loss; 1000],
                ..Default::default()
            })
            .unwrap(),
        ])
        .unwrap();
        let sequence = AASequence::parse_with_registry("AGK(Bulk)GA", &registry).unwrap();
        let generator = TheoreticalSpectrumGenerator {
            add_losses: true,
            ..Default::default()
        };
        let mut work = 50_000_000;
        let mut bytes = 100_000;
        let error = generator
            .generate_for_proforma(&sequence, 1, 1, &mut work, &mut bytes)
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("ProForma theoretical generation resource limit")
        );
        assert!(work < 50_000_000 && bytes < 100_000);
    }
}
