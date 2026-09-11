// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Configurable peptide tandem mass spectra, ported from OpenMS4-core `7c029e8`.
//! Covers terminal ion series, neutral losses, precursors, coarse/fine isotope
//! envelopes, and aligned annotations. See `docs/THEORETICAL_SPECTRA.md`.

use super::fine_isotopes::FineIsotopeWork;
use super::isotopes::CoarseIsotopeWork;
use super::{
    AASequence, CoarseIsotopePatternGenerator, CoarseMassMode, EmpiricalFormula,
    FineIsotopePatternGenerator, FineIsotopeStop, PROTON_MASS_U,
};
use crate::kernel::{DataArray, MSSpectrum, Peak1D, Precursor, SpectrumType};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

#[cfg(test)]
#[path = "theoretical_budget_tests.rs"]
mod budget_tests;
#[path = "theoretical_helpers.rs"]
mod helpers;
#[path = "theoretical_internal.rs"]
mod internal;
#[path = "theoretical_presets.rs"]
mod presets;
#[path = "theoretical_proforma_budget.rs"]
mod proforma_budget;

pub const MAX_THEORETICAL_PEAKS: usize = 100_000;
pub const MAX_THEORETICAL_RESIDUES: usize = 4096;
const MAX_RESIDUE_WORK: usize = 10_000_000;

/// Shared across every candidate and envelope in an annotation operation.
/// Residue visits are precharged conservatively; isotope counters charge the
/// actual convolution products or configuration-search operations.
pub(crate) struct TheoreticalGenerationWork {
    residue_remaining: usize,
    loss_remaining: usize,
    fine: FineIsotopeWork,
    coarse: CoarseIsotopeWork,
}
impl Default for TheoreticalGenerationWork {
    fn default() -> Self {
        Self {
            residue_remaining: 50_000_000,
            loss_remaining: MAX_RESIDUE_WORK,
            fine: FineIsotopeWork::default(),
            coarse: CoarseIsotopeWork::default(),
        }
    }
}
impl TheoreticalGenerationWork {
    fn consume_residues(
        &mut self,
        settings: &TheoreticalSpectrumGenerator,
        residues: usize,
        min_charge: u8,
        max_charge: u8,
    ) -> Result<()> {
        // Same quadratic fragment bound as per-call validation. Internal
        // fragments copy at most ten residues for each of two series/charges.
        let charges = usize::from(max_charge.saturating_sub(min_charge)) + 1;
        let internal = if settings.add_internal_fragments {
            200 * charges
        } else {
            0
        };
        let visits = residues
            .checked_mul(residues)
            .and_then(|n| n.checked_mul(settings.ion_series.len()))
            .and_then(|n| {
                residues
                    .checked_mul(internal + settings.ion_series.len() + charges + 1)
                    .and_then(|extra| n.checked_add(extra))
            })
            .ok_or_else(|| invalid("cumulative theoretical residue work overflow"))?;
        self.residue_remaining = self
            .residue_remaining
            .checked_sub(visits)
            .ok_or_else(|| invalid("cumulative theoretical residue work limit exceeded"))?;
        Ok(())
    }
    fn consume_losses(&mut self, visits: usize) -> Result<()> {
        self.loss_remaining = self
            .loss_remaining
            .checked_sub(visits)
            .ok_or_else(|| invalid("cumulative theoretical loss work limit exceeded"))?;
        Ok(())
    }
}
const ION_NAMES: &str = "IonNames";
const CHARGES: &str = "Charges";

/// Supported series. ZPlusOne/ZPlusTwo are the source's z. and z' variants.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum TheoreticalIonSeries {
    A,
    B,
    C,
    X,
    Y,
    Z,
    ZPlusOne,
    ZPlusTwo,
}

impl TheoreticalIonSeries {
    pub fn label(self) -> &'static str {
        match self {
            Self::A => "a",
            Self::B => "b",
            Self::C => "c",
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
            Self::ZPlusOne => "z.",
            Self::ZPlusTwo => "z'",
        }
    }
    pub fn is_prefix(self) -> bool {
        matches!(self, Self::A | Self::B | Self::C)
    }
    /// Elemental delta relative to the full retained peptide, including its H2O.
    pub fn formula_delta(self) -> EmpiricalFormula {
        let text = match self {
            Self::A => "C-1H-2O-2",
            Self::B => "H-2O-1",
            Self::C => "H1N1O-1",
            Self::X => "C1H-2O1",
            Self::Y => "",
            Self::Z => "H-3N-1",
            Self::ZPlusOne => "H-2N-1",
            Self::ZPlusTwo => "H-1N-1",
        };
        EmpiricalFormula::parse(text).expect("fixed ion conversion formula")
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum TheoreticalIsotopeModel {
    #[default]
    None,
    /// Number of low-mass isotope bins, normalized within each envelope.
    Coarse { max_peaks: usize },
    /// Fine configurations covering 1 minus this probability; no normalization.
    Fine { unexplained_probability: f64 },
}

/// Per-series intact peak intensity and independent precursor/loss intensities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TheoreticalIonIntensities {
    pub a: f32,
    pub b: f32,
    pub c: f32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub precursor: f32,
    pub precursor_water_loss: f32,
    pub precursor_ammonia_loss: f32,
}
impl Default for TheoreticalIonIntensities {
    fn default() -> Self {
        Self {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            x: 1.0,
            y: 1.0,
            z: 1.0,
            precursor: 1.0,
            precursor_water_loss: 1.0,
            precursor_ammonia_loss: 1.0,
        }
    }
}
impl TheoreticalIonIntensities {
    fn for_series(self, series: TheoreticalIonSeries) -> f32 {
        match series {
            TheoreticalIonSeries::A => self.a,
            TheoreticalIonSeries::B => self.b,
            TheoreticalIonSeries::C => self.c,
            TheoreticalIonSeries::X => self.x,
            TheoreticalIonSeries::Y => self.y,
            _ => self.z,
        }
    }
}

/// Theoretical spectrum settings. Defaults match OpenMS's b/y spectrum.
#[derive(Clone, Debug, PartialEq)]
pub struct TheoreticalSpectrumGenerator {
    pub ion_series: Vec<TheoreticalIonSeries>,
    pub intensities: TheoreticalIonIntensities,
    pub isotope_model: TheoreticalIsotopeModel,
    pub add_first_prefix_ion: bool,
    pub add_losses: bool,
    /// Water loss at prefix termini; water and ammonia loss at suffix termini.
    /// Requires add_losses and no isotope envelopes.
    pub add_terminal_losses: bool,
    pub relative_loss_intensity: f32,
    /// Includes intact, water-loss and ammonia-loss precursor peaks, independently
    /// of add_losses. Impossible loss formulas are skipped.
    pub add_precursor_peaks: bool,
    /// Otherwise precursor peaks use max_charge, independently of metadata charge.
    pub add_all_precursor_charges: bool,
    /// Both internal b/a fragments, lengths 2–10, excluding peptide termini.
    /// Independent of terminal ion selection; these peaks have no isotope envelope.
    pub add_internal_fragments: bool,
    /// Fixed singly charged peaks for unmodified P/C/L/H/F/Y/W residues.
    pub add_abundant_immonium_ions: bool,
    pub add_metainfo: bool,
    pub sort_by_position: bool,
}
impl Default for TheoreticalSpectrumGenerator {
    fn default() -> Self {
        Self {
            ion_series: vec![TheoreticalIonSeries::B, TheoreticalIonSeries::Y],
            intensities: TheoreticalIonIntensities::default(),
            isotope_model: TheoreticalIsotopeModel::None,
            add_first_prefix_ion: false,
            add_losses: false,
            add_terminal_losses: false,
            relative_loss_intensity: 0.1,
            add_precursor_peaks: false,
            add_all_precursor_charges: false,
            add_internal_fragments: false,
            add_abundant_immonium_ions: false,
            add_metainfo: false,
            sort_by_position: true,
        }
    }
}

impl TheoreticalSpectrumGenerator {
    /// Generate a centroid MS2 spectrum. None infers precursor charge max_charge+1.
    /// Empty peptides return an empty default spectrum after validating settings.
    pub fn generate(
        &self,
        peptide: &AASequence,
        min_charge: u8,
        max_charge: u8,
        precursor_charge: Option<u16>,
    ) -> Result<MSSpectrum> {
        self.generate_with_work(
            peptide,
            min_charge,
            max_charge,
            precursor_charge,
            &mut TheoreticalGenerationWork::default(),
        )
    }

    pub(crate) fn generate_with_work(
        &self,
        peptide: &AASequence,
        min_charge: u8,
        max_charge: u8,
        precursor_charge: Option<u16>,
        work: &mut TheoreticalGenerationWork,
    ) -> Result<MSSpectrum> {
        if self.ion_series.len() > 8 {
            return Err(invalid("duplicate theoretical ion series"));
        }
        work.consume_residues(self, peptide.len(), min_charge, max_charge)?;
        let precursor_charge = self.validate(peptide, min_charge, max_charge, precursor_charge)?;
        if peptide.is_empty() {
            return Ok(MSSpectrum::default());
        }
        work.consume_losses(self.validate_loss_work(peptide, min_charge, max_charge)?)?;
        let mut loss_entries = MAX_THEORETICAL_PEAKS;
        let mut spectrum = MSSpectrum::default();
        let mut templates = BTreeMap::new();
        // Source insertion order, independent of option-list ordering.
        let ordered = [
            TheoreticalIonSeries::B,
            TheoreticalIonSeries::Y,
            TheoreticalIonSeries::A,
            TheoreticalIonSeries::C,
            TheoreticalIonSeries::X,
            TheoreticalIonSeries::Z,
            TheoreticalIonSeries::ZPlusOne,
            TheoreticalIonSeries::ZPlusTwo,
        ];
        for series in ordered {
            if !self.ion_series.contains(&series) {
                continue;
            }
            let mut ladder = Vec::new();
            let first = if series.is_prefix() && !self.add_first_prefix_ion {
                2
            } else {
                1
            };
            let delta = series.formula_delta();
            for ordinal in first..peptide.len() {
                let fragment = if series.is_prefix() {
                    peptide.prefix(ordinal)?
                } else {
                    peptide.suffix(ordinal)?
                };
                let formula = match fragment.formula() {
                    Ok(formula) => Some(formula.checked_add(&delta)?),
                    Err(Error::Unsupported(_)) if !self.needs_formula() => None,
                    Err(error) => return Err(error),
                };
                if formula.as_ref().is_some_and(|f| !physical_formula(f)) {
                    return Err(invalid("ion series produces a negative atom count"));
                }
                let losses = if self.add_losses {
                    fragment_losses(
                        &fragment,
                        series,
                        self.add_terminal_losses,
                        0,
                        &mut loss_entries,
                    )?
                } else {
                    Vec::new()
                };
                ladder.push(IonTemplate {
                    name: format!("{}{ordinal}", series.label()),
                    mono_mass: fragment.mono_mass()? + delta.mono_mass(),
                    formula,
                    losses,
                });
            }
            templates.insert(series, ladder);
        }
        let envelope_peaks = match self.isotope_model {
            TheoreticalIsotopeModel::None => 1,
            TheoreticalIsotopeModel::Coarse { max_peaks } => max_peaks,
            // Fine envelope lengths are discovered during bounded enumeration.
            // This reserves their minimum; push_peak checks the actual total.
            TheoreticalIsotopeModel::Fine {
                unexplained_probability,
            } => usize::from(unexplained_probability < 1.0),
        };
        let charge_count = usize::from(max_charge - min_charge) + 1;
        let mut count = 0_usize;
        for template in templates.values().flatten() {
            count = count
                .checked_add(1 + template.losses.len())
                .ok_or_else(|| invalid("theoretical peak count overflow"))?;
        }
        count = count
            .checked_mul(charge_count)
            .ok_or_else(|| invalid("theoretical peak count overflow"))?;
        if self.add_precursor_peaks {
            count = count
                .checked_add(
                    3 * if self.add_all_precursor_charges {
                        charge_count
                    } else {
                        1
                    },
                )
                .ok_or_else(|| invalid("theoretical peak count overflow"))?;
        }
        count = count
            .checked_mul(envelope_peaks)
            .ok_or_else(|| invalid("theoretical isotope count overflow"))?;
        // Internal and immonium peaks never expand into isotope envelopes.
        // Internal losses are bounded again by push_peak as they are discovered.
        if self.add_internal_fragments {
            let internal_count: usize = (1..peptide.len().saturating_sub(3))
                .map(|start| (peptide.len() - 1 - start).min(10) - 1)
                .sum();
            count = count
                .checked_add(internal_count * 2 * charge_count)
                .ok_or_else(|| invalid("theoretical internal count overflow"))?;
        }
        if self.add_abundant_immonium_ions {
            count = count
                .checked_add(7)
                .ok_or_else(|| invalid("theoretical immonium count overflow"))?;
        }
        if count > MAX_THEORETICAL_PEAKS {
            return Err(invalid("requested theoretical spectrum exceeds peak limit"));
        }
        let mut rows = Vec::with_capacity(count);
        let mut emitter = Emitter::new(self, &mut work.fine, &mut work.coarse)?;
        for charge in min_charge..=max_charge {
            for series in ordered {
                let Some(ladder) = templates.get(&series) else {
                    continue;
                };
                let intensity = self.intensities.for_series(series);
                for ion in ladder {
                    let annotation = if self.isotope_model == TheoreticalIsotopeModel::None {
                        format!("{}{}", ion.name, "+".repeat(usize::from(charge)))
                    } else {
                        ion.name.clone()
                    };
                    emitter.emit(
                        &mut rows,
                        ion.formula.as_ref(),
                        ion.mono_mass,
                        charge,
                        f64::from(intensity),
                        &annotation,
                    )?;
                }
                for ion in ladder {
                    for loss in &ion.losses {
                        let loss_formula = ion
                            .formula
                            .as_ref()
                            .map(|formula| formula.checked_sub(loss))
                            .transpose()?;
                        if loss_formula.as_ref().is_some_and(|f| !physical_formula(f)) {
                            continue;
                        }
                        let annotation =
                            format!("{}-{}{}", ion.name, loss, "+".repeat(usize::from(charge)));
                        emitter.emit(
                            &mut rows,
                            loss_formula.as_ref(),
                            ion.mono_mass - loss.mono_mass(),
                            charge,
                            f64::from(intensity) * f64::from(self.relative_loss_intensity),
                            &annotation,
                        )?;
                    }
                }
            }
        }
        if self.add_internal_fragments {
            self.internal_fragment_peaks(
                &mut rows,
                peptide,
                min_charge,
                max_charge,
                &mut loss_entries,
            )?;
        }
        if self.add_precursor_peaks {
            let start = if self.add_all_precursor_charges {
                min_charge
            } else {
                max_charge
            };
            for charge in start..=max_charge {
                self.precursor_peaks(&mut emitter, &mut rows, peptide, charge)?;
            }
        }
        if self.add_abundant_immonium_ions {
            self.emit_immonium(&mut rows, peptide)?;
        }
        if self.sort_by_position {
            rows.sort_by(|a, b| a.peak.mz.total_cmp(&b.peak.mz));
        }
        spectrum.peaks = rows.iter().map(|r| r.peak).collect();
        if self.add_metainfo {
            spectrum.string_data_arrays.push(DataArray::new(
                ION_NAMES,
                rows.iter().map(|r| r.name.clone()).collect(),
            ));
            spectrum.integer_data_arrays.push(DataArray::new(
                CHARGES,
                rows.iter().map(|r| i32::from(r.charge)).collect(),
            ));
        }
        spectrum.ms_level = 2;
        spectrum.spectrum_type = SpectrumType::Centroid;
        spectrum.precursors.push(Precursor::new(
            peptide.mz(i32::from(precursor_charge))?,
            i32::from(precursor_charge),
        ));
        spectrum.validate()?;
        Ok(spectrum)
    }

    /// Append atomically, keeping existing peaks, precursor entries and metadata.
    /// Named IonNames/Charges arrays are extended (existing unannotated peaks get
    /// empty names/charge zero). Unrelated populated per-peak arrays are rejected:
    /// there is no defined value for their new theoretical peaks. Empty placeholders
    /// are retained. Any error leaves the complete target unchanged.
    pub fn append_to(
        &self,
        spectrum: &mut MSSpectrum,
        peptide: &AASequence,
        min_charge: u8,
        max_charge: u8,
        precursor_charge: Option<u16>,
    ) -> Result<()> {
        // Acquisition settings are owned by the appended record; meter their
        // copy separately from the unchanged theoretical-generation work limit.
        spectrum.acquisition_with_budget(&mut 50_000_000, &mut (256 * 1024 * 1024))?;
        spectrum.validate()?;
        let names = unique_named_array(&spectrum.string_data_arrays, ION_NAMES)?;
        let charges = unique_named_array(&spectrum.integer_data_arrays, CHARGES)?;
        let mut settings = self.clone();
        if names.is_some() || charges.is_some() {
            settings.add_metainfo = true;
        }
        let addition = settings.generate(peptide, min_charge, max_charge, precursor_charge)?;
        if peptide.is_empty() {
            return Ok(());
        }
        if spectrum
            .len()
            .checked_add(addition.len())
            .is_none_or(|n| n > MAX_THEORETICAL_PEAKS)
        {
            return Err(invalid("combined theoretical spectrum exceeds peak limit"));
        }
        if !addition.is_empty()
            && (spectrum
                .float_data_arrays
                .iter()
                .any(|a| !a.data.is_empty())
                || spectrum
                    .string_data_arrays
                    .iter()
                    .any(|a| a.name != ION_NAMES && !a.data.is_empty())
                || spectrum
                    .integer_data_arrays
                    .iter()
                    .any(|a| a.name != CHARGES && !a.data.is_empty()))
        {
            return Err(Error::Unsupported(
                "cannot append peaks to unrelated populated annotation arrays".into(),
            ));
        }
        let mut combined = spectrum.clone();
        let old_length = combined.len();
        combined.peaks.extend(addition.peaks);
        if settings.add_metainfo {
            let new_names = addition
                .string_data_arrays
                .into_iter()
                .next()
                .expect("requested annotations")
                .data;
            let new_charges = addition
                .integer_data_arrays
                .into_iter()
                .next()
                .expect("requested annotations")
                .data;
            append_annotations(
                &mut combined.string_data_arrays,
                names,
                ION_NAMES,
                old_length,
                new_names,
            );
            append_annotations(
                &mut combined.integer_data_arrays,
                charges,
                CHARGES,
                old_length,
                new_charges,
            );
        }
        combined.precursors.extend(addition.precursors);
        combined.ms_level = 2;
        combined.spectrum_type = SpectrumType::Centroid;
        if settings.sort_by_position {
            combined.sort_by_position()?;
        }
        combined.validate()?;
        *spectrum = combined;
        Ok(())
    }

    fn validate(
        &self,
        peptide: &AASequence,
        min_charge: u8,
        max_charge: u8,
        precursor_charge: Option<u16>,
    ) -> Result<u16> {
        if min_charge == 0 || min_charge > max_charge {
            return Err(invalid("fragment charges require 1 <= min <= max"));
        }
        let precursor_charge = precursor_charge.unwrap_or(u16::from(max_charge) + 1);
        if precursor_charge < u16::from(max_charge) {
            return Err(invalid(
                "precursor charge must be at least the maximum fragment charge",
            ));
        }
        let unique: BTreeSet<_> = self.ion_series.iter().copied().collect();
        if unique.len() != self.ion_series.len() {
            return Err(invalid("duplicate theoretical ion series"));
        }
        let n = peptide.len();
        if n > MAX_THEORETICAL_RESIDUES
            || n.checked_mul(n)
                .and_then(|work| work.checked_mul(self.ion_series.len()))
                .is_none_or(|work| work > MAX_RESIDUE_WORK)
        {
            return Err(invalid("peptide exceeds theoretical fragment work limit"));
        }
        if n == 1
            && self
                .ion_series
                .iter()
                .any(|s| matches!(s, TheoreticalIonSeries::C | TheoreticalIonSeries::X))
        {
            return Err(invalid("c/x ion generation requires at least two residues"));
        }
        let i = self.intensities;
        for value in [
            i.a,
            i.b,
            i.c,
            i.x,
            i.y,
            i.z,
            i.precursor,
            i.precursor_water_loss,
            i.precursor_ammonia_loss,
            self.relative_loss_intensity,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(invalid(
                    "theoretical intensities must be finite and nonnegative",
                ));
            }
        }
        if self.relative_loss_intensity > 1.0 {
            return Err(invalid("relative loss intensity cannot exceed one"));
        }
        if self.add_terminal_losses
            && (!self.add_losses || self.isotope_model != TheoreticalIsotopeModel::None)
        {
            return Err(Error::Unsupported(
                "terminal-loss additions require add_losses without isotope envelopes".into(),
            ));
        }
        if let TheoreticalIsotopeModel::Coarse { max_peaks } = self.isotope_model {
            if max_peaks == 0 || max_peaks > MAX_THEORETICAL_PEAKS {
                return Err(invalid("invalid theoretical isotope peak count"));
            }
        }
        if let TheoreticalIsotopeModel::Fine {
            unexplained_probability,
        } = self.isotope_model
        {
            if !unexplained_probability.is_finite()
                || !(0.0..=1.0).contains(&unexplained_probability)
            {
                return Err(invalid(
                    "fine isotope unexplained probability must be in [0, 1]",
                ));
            }
        }
        if !peptide.is_empty() {
            peptide.mono_mass()?;
            if self.needs_formula() {
                peptide.formula()?;
            }
        }
        Ok(precursor_charge)
    }

    // Monoisotopic intact ions can use an observed mass without assigning atoms.
    // Isotope probabilities and the source's formula-based precursor losses
    // need a known composition. Fragment losses use declared residue losses.
    fn needs_formula(&self) -> bool {
        self.add_precursor_peaks
            || (self.isotope_model != TheoreticalIsotopeModel::None && !self.ion_series.is_empty())
    }

    // Charge all declarations before scanning or deduplicating them. Repeated
    // custom losses can have small output but otherwise arbitrarily large work.
    fn validate_loss_work(
        &self,
        peptide: &AASequence,
        min_charge: u8,
        max_charge: u8,
    ) -> Result<usize> {
        if !self.add_losses {
            return Ok(0);
        }
        let costs = (0..peptide.len())
            .map(|index| {
                let declarations = peptide
                    .residue_modification(index)?
                    .map_or(3, |m| m.known().map_or(0, |m| m.neutral_losses().len()));
                declarations
                    .checked_add(1)
                    .ok_or_else(|| invalid("neutral loss work overflow"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut total = 0_usize;
        let mut add = |cost: usize, visits: usize| -> Result<()> {
            total = cost
                .checked_mul(visits)
                .and_then(|n| total.checked_add(n))
                .filter(|&n| n <= MAX_RESIDUE_WORK)
                .ok_or_else(|| invalid("neutral loss declaration work exceeds limit"))?;
            Ok(())
        };
        for series in &self.ion_series {
            let first = if series.is_prefix() && !self.add_first_prefix_ion {
                2
            } else {
                1
            };
            for (index, &cost) in costs.iter().enumerate() {
                let first_included = if series.is_prefix() {
                    index + 1
                } else {
                    peptide.len() - index
                };
                add(
                    cost,
                    peptide.len().saturating_sub(first.max(first_included)),
                )?;
            }
        }
        if self.add_internal_fragments {
            let multiplier = 2 * (usize::from(max_charge - min_charge) + 1);
            for start in 1..peptide.len().saturating_sub(3) {
                let end = (start + 10).min(peptide.len() - 1);
                for (index, &cost) in costs.iter().enumerate().take(end).skip(start + 1) {
                    add(cost, (end - index) * multiplier)?;
                }
            }
        }
        Ok(total)
    }

    fn precursor_peaks(
        &self,
        emitter: &mut Emitter<'_>,
        rows: &mut Vec<PeakRow>,
        peptide: &AASequence,
        charge: u8,
    ) -> Result<()> {
        let suffix = "+".repeat(usize::from(charge));
        let adduct = if charge == 1 {
            "H".into()
        } else {
            format!("{charge}H")
        };
        let formula = peptide.formula()?;
        emitter.emit(
            rows,
            Some(&formula),
            peptide.mono_mass()?,
            charge,
            f64::from(self.intensities.precursor),
            &format!("[M+{adduct}]{suffix}"),
        )?;
        for (loss, label, intensity) in [
            ("H2O", "H2O", self.intensities.precursor_water_loss),
            ("NH3", "NH3", self.intensities.precursor_ammonia_loss),
        ] {
            let loss_formula = formula.checked_sub(&EmpiricalFormula::parse(loss)?)?;
            if physical_formula(&loss_formula) {
                // The source uses formula-derived masses for precursor losses,
                // even though its intact precursor retains tabulated terminal deltas.
                emitter.emit(
                    rows,
                    Some(&loss_formula),
                    loss_formula.mono_mass(),
                    charge,
                    f64::from(intensity),
                    &format!("[M+{adduct}-{label}]{suffix}"),
                )?;
            }
        }
        Ok(())
    }
}

struct IonTemplate {
    name: String,
    formula: Option<EmpiricalFormula>,
    mono_mass: f64,
    losses: Vec<EmpiricalFormula>,
}
struct PeakRow {
    peak: Peak1D,
    name: String,
    charge: u8,
}
struct Emitter<'a> {
    isotope_generator: Option<EnvelopeGenerator>,
    fine_work: &'a mut FineIsotopeWork,
    coarse_work: &'a mut CoarseIsotopeWork,
}
enum EnvelopeGenerator {
    Coarse(CoarseIsotopePatternGenerator),
    Fine(FineIsotopePatternGenerator),
}
impl<'a> Emitter<'a> {
    fn new(
        settings: &TheoreticalSpectrumGenerator,
        fine_work: &'a mut FineIsotopeWork,
        coarse_work: &'a mut CoarseIsotopeWork,
    ) -> Result<Self> {
        let isotope_generator = match settings.isotope_model {
            TheoreticalIsotopeModel::None => None,
            TheoreticalIsotopeModel::Coarse { max_peaks } => Some(EnvelopeGenerator::Coarse(
                CoarseIsotopePatternGenerator::new(Some(max_peaks), CoarseMassMode::Approximate)?,
            )),
            TheoreticalIsotopeModel::Fine {
                unexplained_probability,
            } => Some(EnvelopeGenerator::Fine(FineIsotopePatternGenerator {
                stop: FineIsotopeStop::UnexplainedProbability(unexplained_probability),
            })),
        };
        Ok(Self {
            isotope_generator,
            fine_work,
            coarse_work,
        })
    }
    fn emit(
        &mut self,
        rows: &mut Vec<PeakRow>,
        formula: Option<&EmpiricalFormula>,
        mono_mass: f64,
        charge: u8,
        intensity: f64,
        annotation: &str,
    ) -> Result<()> {
        if let Some(generator) = &self.isotope_generator {
            let formula = formula
                .ok_or_else(|| invalid("isotope generation requires an elemental formula"))?;
            // Match TSG's explicit neutral H adduct, not the generator's implicit
            // protonation behavior. It retains the electron mass in the envelope.
            let adduct = EmpiricalFormula::parse("H")?.checked_scale(i32::from(charge))?;
            let ion = formula.checked_add(&adduct)?.with_charge(0);
            let envelope = match generator {
                EnvelopeGenerator::Coarse(generator) => {
                    generator.run_with_work(&ion, self.coarse_work)?
                }
                EnvelopeGenerator::Fine(generator) => {
                    generator.run_with_work(&ion, self.fine_work)?
                }
            };
            for peak in envelope.peaks() {
                push_peak(
                    rows,
                    peak.mass / f64::from(charge),
                    (intensity * peak.probability) as f32,
                    annotation,
                    charge,
                )?;
            }
        } else {
            push_peak(
                rows,
                (mono_mass + f64::from(charge) * PROTON_MASS_U) / f64::from(charge),
                intensity as f32,
                annotation,
                charge,
            )?;
        }
        Ok(())
    }
}

fn push_peak(
    rows: &mut Vec<PeakRow>,
    mz: f64,
    intensity: f32,
    name: &str,
    charge: u8,
) -> Result<()> {
    if rows.len() >= MAX_THEORETICAL_PEAKS {
        return Err(invalid("theoretical spectrum exceeds peak limit"));
    }
    if !mz.is_finite() || mz < 0.0 || !intensity.is_finite() || intensity < 0.0 {
        return Err(invalid(
            "theoretical peak must have finite, nonnegative m/z and intensity",
        ));
    }
    rows.push(PeakRow {
        peak: Peak1D::new(mz, intensity),
        name: name.into(),
        charge,
    });
    Ok(())
}

fn fragment_losses(
    fragment: &AASequence,
    series: TheoreticalIonSeries,
    terminal: bool,
    skip: usize,
    remaining_entries: &mut usize,
) -> Result<Vec<EmpiricalFormula>> {
    let mut losses = BTreeMap::new();
    let mut add = |formula: &EmpiricalFormula| -> Result<()> {
        if let std::collections::btree_map::Entry::Vacant(entry) = losses.entry(formula.to_string())
        {
            *remaining_entries = remaining_entries
                .checked_sub(1)
                .ok_or_else(|| invalid("neutral loss template storage exceeds limit"))?;
            entry.insert(formula.clone());
        }
        Ok(())
    };
    for (index, residue) in fragment.as_str().bytes().enumerate().skip(skip) {
        if let Some(modification) = fragment.residue_modification(index)? {
            // OpenMS replaces the unmodified residue's losses when modified.
            for loss in modification
                .known()
                .into_iter()
                .flat_map(|m| m.neutral_losses())
            {
                if !loss.formula().is_empty() {
                    add(loss.formula())?;
                }
            }
        } else {
            let formulas: &[&str] = match residue {
                b'D' | b'E' | b'S' | b'T' => &["H2O"],
                b'K' | b'N' | b'Q' => &["NH3"],
                b'R' => &["NH3", "C1H2N2", "C1H2N1O1"],
                _ => &[],
            };
            for text in formulas {
                let formula = EmpiricalFormula::parse(text)?;
                add(&formula)?;
            }
        }
    }
    if terminal {
        let water = EmpiricalFormula::parse("H2O")?;
        add(&water)?;
        if !series.is_prefix() {
            let ammonia = EmpiricalFormula::parse("NH3")?;
            add(&ammonia)?;
        }
    }
    Ok(losses.into_values().collect())
}

fn physical_formula(formula: &EmpiricalFormula) -> bool {
    formula.atoms.values().all(|&count| count >= 0)
}
fn unique_named_array<T>(arrays: &[DataArray<T>], name: &str) -> Result<Option<usize>> {
    let indices: Vec<_> = arrays
        .iter()
        .enumerate()
        .filter(|(_, a)| a.name == name)
        .map(|(i, _)| i)
        .collect();
    if indices.len() > 1 {
        return Err(invalid(format!("duplicate {name} annotation arrays")));
    }
    Ok(indices.first().copied())
}
fn append_annotations<T: Default + Clone>(
    arrays: &mut Vec<DataArray<T>>,
    index: Option<usize>,
    name: &str,
    old_length: usize,
    addition: Vec<T>,
) {
    if let Some(index) = index {
        if arrays[index].data.is_empty() {
            arrays[index].data.resize(old_length, T::default());
        }
        arrays[index].data.extend(addition);
    } else {
        let mut data = vec![T::default(); old_length];
        data.extend(addition);
        arrays.push(DataArray::new(name, data));
    }
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
