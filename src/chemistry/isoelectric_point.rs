// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Peptide charge and pI using the pinned OpenMS proteomics pKa scales.
//! Parent residue codes determine sidechain charge. Any terminal annotation
//! suppresses that terminal group; PTM-specific pKa shifts are not modeled.

use super::AASequence;
use crate::{Error, Result};

/// Maximum peptide length inspected by either calculation.
pub const MAX_PI_RESIDUES: usize = 1_000_000;
/// Shared budget of residue/terminal visits across all charge evaluations.
pub const MAX_PI_WORK: usize = 50_000_000;
/// Additional guard against impractically fine bisection requests.
pub const MAX_PI_ITERATIONS: usize = 128;

/// Published pKa tables used by the source IsoelectricPoint utility.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ProteomicsPkaScale {
    #[default]
    Lehninger,
    Emboss,
    Sillero,
    Bjellqvist,
}

/// Charge-model options. The tolerance is used only by [`Self::compute_pi`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsoelectricPoint {
    pub scale: ProteomicsPkaScale,
    /// Positive finite endpoint-charge and pH-interval tolerance (source default 1e-4).
    pub tolerance: f64,
}

impl Default for IsoelectricPoint {
    fn default() -> Self {
        Self {
            scale: ProteomicsPkaScale::default(),
            tolerance: 1e-4,
        }
    }
}

impl IsoelectricPoint {
    /// Net charge at any finite pH, including outside the pI search interval.
    /// B/Z/X/J have neutral sidechains, U uses acidic pKa 5.73, and O errors.
    /// Empty/oversized peptides and nonfinite pH error. Unused tolerance is not
    /// validated. Infinite powers giving finite saturated charges are allowed.
    pub fn compute_charge(&self, sequence: &AASequence, ph: f64) -> Result<f64> {
        let mut remaining = MAX_PI_WORK;
        self.charge_with_work(sequence, ph, &mut remaining)
    }

    /// Bisect `[0, 14]` with source endpoint priority and interval-width stopping.
    /// If charge does not cross zero, return the boundary with smaller absolute
    /// charge (ties choose zero). Any endpoint within tolerance returns first,
    /// checking zero before fourteen. An exactly zero midpoint advances high.
    /// Invalid tolerance, exhausted shared work/iterations or a stalled floating
    /// midpoint returns an error instead of a partial estimate or endless loop.
    pub fn compute_pi(&self, sequence: &AASequence) -> Result<f64> {
        let mut remaining = MAX_PI_WORK;
        self.pi_with_work(sequence, &mut remaining)
    }

    fn pi_with_work(&self, sequence: &AASequence, remaining: &mut usize) -> Result<f64> {
        validate_sequence(sequence)?;
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err(invalid("pI tolerance must be finite and positive"));
        }
        let mut low = 0.0;
        let mut high = 14.0;
        let charge_low = self.charge_with_work(sequence, low, remaining)?;
        if charge_low.abs() <= self.tolerance {
            return Ok(low);
        }
        let charge_high = self.charge_with_work(sequence, high, remaining)?;
        if charge_high.abs() <= self.tolerance {
            return Ok(high);
        }
        if (charge_low > 0.0 && charge_high > 0.0) || (charge_low < 0.0 && charge_high < 0.0) {
            return Ok(if charge_low.abs() <= charge_high.abs() {
                low
            } else {
                high
            });
        }

        let mut iterations = 0;
        while high - low > self.tolerance {
            if iterations >= MAX_PI_ITERATIONS {
                return Err(invalid("pI iteration limit exceeded"));
            }
            let midpoint = (low + high) / 2.0;
            if midpoint == low || midpoint == high {
                return Err(invalid(
                    "pI interval cannot progress at requested tolerance",
                ));
            }
            let charge = self.charge_with_work(sequence, midpoint, remaining)?;
            if charge > 0.0 {
                low = midpoint;
            } else {
                high = midpoint;
            }
            iterations += 1;
        }
        Ok((low + high) / 2.0)
    }

    fn charge_with_work(
        &self,
        sequence: &AASequence,
        ph: f64,
        remaining: &mut usize,
    ) -> Result<f64> {
        validate_sequence(sequence)?;
        if !ph.is_finite() {
            return Err(invalid("charge pH must be finite"));
        }
        // Charge even neutral residues and omitted termini before traversing.
        *remaining = remaining
            .checked_sub(sequence.len() + 2)
            .ok_or_else(|| invalid("pI work limit exceeded"))?;
        let [mut nterm, mut cterm, d, e, c, y, h, k, r] = match self.scale {
            ProteomicsPkaScale::Lehninger => {
                [9.69, 2.34, 3.65, 4.25, 8.18, 10.07, 6.00, 10.53, 12.48]
            }
            ProteomicsPkaScale::Emboss => [8.6, 3.6, 3.9, 4.1, 8.5, 10.1, 6.5, 10.8, 12.5],
            ProteomicsPkaScale::Sillero => [8.2, 3.2, 4.0, 4.5, 9.0, 10.0, 6.4, 10.4, 12.0],
            ProteomicsPkaScale::Bjellqvist => [7.50, 3.55, 4.05, 4.45, 9.0, 10.0, 5.98, 10.0, 12.0],
        };
        let residues = sequence.as_str().as_bytes();
        if self.scale == ProteomicsPkaScale::Bjellqvist {
            nterm = match residues[0] {
                b'A' => 7.59,
                b'R' | b'E' => 7.70,
                b'C' => 8.00,
                b'M' => 7.00,
                b'P' => 8.36,
                b'S' => 6.93,
                b'T' => 6.82,
                b'V' => 7.44,
                // Other canonical residues and unrecognized codes use 7.50.
                _ => nterm,
            };
            cterm = match residues[residues.len() - 1] {
                b'D' => 4.55,
                b'E' => 4.75,
                _ => cterm,
            };
        }
        let mut charge = 0.0;
        if sequence.n_terminal_modification().is_none() {
            charge += basic(ph, nterm);
        }
        if sequence.c_terminal_modification().is_none() {
            charge += acidic(ph, cterm);
        }
        // Source adds each contribution in residue order, without count pooling.
        for residue in residues {
            match residue {
                b'D' => charge += acidic(ph, d),
                b'E' => charge += acidic(ph, e),
                b'C' => charge += acidic(ph, c),
                b'Y' => charge += acidic(ph, y),
                b'U' => charge += acidic(ph, 5.73),
                b'H' => charge += basic(ph, h),
                b'K' => charge += basic(ph, k),
                b'R' => charge += basic(ph, r),
                b'O' => {
                    return Err(Error::Unsupported(
                        "pyrrolysine is not supported for pI/charge calculation".into(),
                    ));
                }
                _ => {}
            }
        }
        if !charge.is_finite() {
            return Err(invalid("peptide charge is nonfinite"));
        }
        Ok(charge)
    }
}

fn acidic(ph: f64, pka: f64) -> f64 {
    -1.0 / (1.0 + 10.0_f64.powf(pka - ph))
}
fn basic(ph: f64, pka: f64) -> f64 {
    1.0 / (1.0 + 10.0_f64.powf(ph - pka))
}
fn validate_sequence(sequence: &AASequence) -> Result<()> {
    if sequence.is_empty() {
        return Err(invalid("cannot compute pI/charge of an empty sequence"));
    }
    if sequence.len() > MAX_PI_RESIDUES {
        return Err(invalid("pI residue limit exceeded"));
    }
    Ok(())
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn endpoint_and_midpoint_evaluations_share_one_precharged_budget() {
        let peptide = AASequence::parse("A").unwrap();
        let model = IsoelectricPoint::default();
        // Three visits each: both endpoints and first midpoint consume nine.
        // The next evaluation must fail rather than restart a fresh budget.
        let mut remaining = 9;
        let error = model.pi_with_work(&peptide, &mut remaining).unwrap_err();
        assert!(error.to_string().contains("pI work limit"));
        assert_eq!(remaining, 0);

        let mut remaining = 2;
        assert!(
            model
                .charge_with_work(&peptide, 7.0, &mut remaining)
                .is_err()
        );
        assert_eq!(remaining, 2); // Failed precharge does not traverse the input.
    }
}
