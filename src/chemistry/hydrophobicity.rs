// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Residue hydrophobicity, GRAVY, rolling profiles and hydrophobic moments.
//!
//! Tables and arithmetic follow OpenMS4-core revision `7c029e8`,
//! [`Residue.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8/src/openms/source/CHEMISTRY/Residue.cpp)
//! and `HydrophobicityProfile.cpp`. Values refer to parent amino acids; sequence
//! annotations do not alter any of these calculations.

use super::AASequence;
use super::aa_index::{canonical_index, check_property_limits, finite, invalid};
use crate::Result;

/// Seven source hydrophobicity scales, without reorientation or normalization.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum HydrophobicityScale {
    /// Kyte–Doolittle (1982); also used by GRAVY.
    #[default]
    KyteDoolittle,
    /// Normalized Eisenberg (1984); also used by hydrophobic moments.
    Eisenberg,
    /// Hopp–Woods (1981).
    HoppWoods,
    /// Bull–Breese (1974).
    BullBreese,
    /// Black–Mould (1991).
    BlackMould,
    /// Guy (1985).
    Guy,
    /// Eisenberg consensus, distinct from normalized Eisenberg.
    EisenbergConsensus,
}

impl HydrophobicityScale {
    /// Every scale in the source enum/table order.
    pub const ALL: [Self; 7] = [
        Self::KyteDoolittle,
        Self::Eisenberg,
        Self::HoppWoods,
        Self::BullBreese,
        Self::BlackMould,
        Self::Guy,
        Self::EisenbergConsensus,
    ];

    /// One uppercase canonical residue's value. The six source sentinel codes
    /// B/J/O/U/X/Z, lowercase characters and all other codes return an error.
    pub fn value(self, residue: char) -> Result<f64> {
        Ok(SCALE_VALUES[self as usize][canonical_index(residue)?])
    }
}

/// Source peptide hydrophobicity calculations with checked owned outputs.
///
/// All functions ignore annotations and require a nonempty canonical parent
/// sequence. Input/output lengths are limited to 1,000,000 residues/values and
/// total charged work to 50,000,000 units before allocation or traversal.
#[derive(Clone, Copy, Debug, Default)]
pub struct HydrophobicityProfile;

impl HydrophobicityProfile {
    /// Sequential mean Kyte–Doolittle hydrophobicity (GRAVY).
    pub fn compute_gravy(sequence: &AASequence) -> Result<f64> {
        let residues = prepare(sequence, sequence.len().checked_mul(2))?;
        let mut sum = 0.0;
        for &residue in residues {
            sum += HydrophobicityScale::KyteDoolittle.value(char::from(residue))?;
        }
        finite(sum / residues.len() as f64)
    }

    /// One value per parent residue. The source default scale is Kyte–Doolittle.
    pub fn compute_profile(sequence: &AASequence, scale: HydrophobicityScale) -> Result<Vec<f64>> {
        let residues = prepare(sequence, sequence.len().checked_mul(2))?;
        let mut output = Vec::with_capacity(residues.len());
        for &residue in residues {
            output.push(scale.value(char::from(residue))?);
        }
        Ok(output)
    }

    /// Rolling means, removing the old value before adding the new one.
    ///
    /// The source defaults are window 7 and Kyte–Doolittle. A positive window
    /// larger than the sequence is clamped to its length; zero is an error.
    pub fn compute_windowed_profile(
        sequence: &AASequence,
        window: usize,
        scale: HydrophobicityScale,
    ) -> Result<Vec<f64>> {
        let width = window_width(sequence, window)?;
        // Full validation plus <=two value lookups per residue.
        let residues = prepare(sequence, sequence.len().checked_mul(3))?;
        let mut output = Vec::with_capacity(residues.len() - width + 1);
        let mut sum = 0.0;
        for &residue in &residues[..width] {
            sum += scale.value(char::from(residue))?;
        }
        output.push(finite(sum / width as f64)?);
        for end in width..residues.len() {
            sum -= scale.value(char::from(residues[end - width]))?;
            sum += scale.value(char::from(residues[end]))?;
            output.push(finite(sum / width as f64)?);
        }
        Ok(output)
    }

    /// Eisenberg normalized hydrophobic moments with window phase reset to zero.
    ///
    /// The source defaults are window 11 and 100 degrees (alpha helices); 160
    /// degrees is the common beta-sheet choice. Every finite angle is accepted
    /// if its source-ordered radian and phase products remain finite. Angles are
    /// not reduced modulo 360. Positive windows are clamped to sequence length.
    ///
    /// Work is precharged as `2*n + 3*w + 2*w*(n-w+1) + (n-w+1)`: validation,
    /// values, phases/trigonometry, two sums and final results. Values and phase
    /// sine/cosine pairs are cached, preserving the source multiplication and
    /// addition order within every window. No residue/subsequence is cloned.
    pub fn compute_hydrophobic_moment(
        sequence: &AASequence,
        window: usize,
        angle_degrees: f64,
    ) -> Result<Vec<f64>> {
        let width = window_width(sequence, window)?;
        if !angle_degrees.is_finite() {
            return Err(invalid("hydrophobic-moment angle must be finite"));
        }
        let outputs = sequence.len() - width + 1;
        let work = sequence
            .len()
            .checked_mul(2)
            .and_then(|work| {
                width
                    .checked_mul(3)
                    .and_then(|setup| work.checked_add(setup))
            })
            .and_then(|work| {
                width
                    .checked_mul(outputs)
                    .and_then(|samples| samples.checked_mul(2))
                    .and_then(|samples| work.checked_add(samples))
            })
            .and_then(|work| work.checked_add(outputs));
        let residues = prepare(sequence, work)?;
        // Constants::PI has exactly the same binary64 value as Rust's PI. Keep
        // multiplication before division: dividing degrees first changes bits.
        let angle_rad = angle_degrees * std::f64::consts::PI / 180.0;
        if !angle_rad.is_finite() || !(angle_rad * (width - 1) as f64).is_finite() {
            return Err(invalid(
                "hydrophobic-moment angle or phase product overflowed",
            ));
        }
        let mut values = Vec::with_capacity(residues.len());
        for &residue in residues {
            values.push(HydrophobicityScale::Eisenberg.value(char::from(residue))?);
        }
        let mut phases = Vec::with_capacity(width);
        for position in 0..width {
            let phase = angle_rad * position as f64;
            phases.push((phase.sin(), phase.cos()));
        }
        let mut output = Vec::with_capacity(outputs);
        for values in values.windows(width) {
            let mut sum_sin = 0.0;
            let mut sum_cos = 0.0;
            for (&hydrophobicity, &(sin, cos)) in values.iter().zip(&phases) {
                sum_sin += hydrophobicity * sin;
                sum_cos += hydrophobicity * cos;
            }
            // Preserve the source square/add/sqrt expression, not hypot.
            sum_sin *= sum_sin;
            sum_cos *= sum_cos;
            output.push(finite((sum_sin + sum_cos).sqrt() / width as f64)?);
        }
        Ok(output)
    }
}

fn window_width(sequence: &AASequence, window: usize) -> Result<usize> {
    if sequence.is_empty() {
        return Err(invalid("hydrophobicity requires a nonempty sequence"));
    }
    if window == 0 {
        return Err(invalid("hydrophobicity window must be positive"));
    }
    Ok(window.min(sequence.len()))
}

fn prepare(sequence: &AASequence, work: Option<usize>) -> Result<&[u8]> {
    check_property_limits(sequence.len(), work)?;
    if sequence.is_empty() {
        return Err(invalid("hydrophobicity requires a nonempty sequence"));
    }
    let residues = sequence.as_str().as_bytes();
    for &residue in residues {
        canonical_index(char::from(residue))?;
    }
    Ok(residues)
}

// Exact Residue.cpp scale values for ACDEFGHIKLMNPQRSTVWY. The source's 999
// placeholders for B/J/O/U/X/Z become checked errors rather than numeric values.
const SCALE_VALUES: [[f64; 20]; 7] = [
    [
        1.800, 2.500, -3.50, -3.50, 2.800, -0.40, -3.20, 4.500, -3.90, 3.800, 1.900, -3.50, -1.60,
        -3.50, -4.50, -0.80, -0.70, 4.200, -0.90, -1.30,
    ],
    [
        0.620, 0.290, -0.90, -0.74, 1.190, 0.480, -0.40, 1.380, -1.50, 1.060, 0.640, -0.780, 0.120,
        -0.85, -2.53, -0.18, -0.05, 1.080, 0.810, 0.260,
    ],
    [
        -0.50, -1.00, 3.000, 3.000, -2.50, 0.000, -0.50, -1.80, 3.000, -1.80, -1.30, 0.200, 0.000,
        0.200, 3.000, 0.300, -0.40, -1.50, -3.40, -2.30,
    ],
    [
        0.610, 0.360, 0.610, 0.510, -1.52, 0.810, 0.690, -1.45, 0.460, -1.65, -0.66, 0.890, -0.17,
        0.970, 0.690, 0.420, 0.290, -0.75, -1.20, -1.43,
    ],
    [
        0.616, 0.680, 0.028, 0.043, 1.000, 0.501, 0.165, 0.943, 0.283, 0.943, 0.738, 0.236, 0.711,
        0.251, 0.000, 0.359, 0.450, 0.825, 0.878, 0.880,
    ],
    [
        0.100, -1.42, 0.780, 0.830, -2.12, 0.330, -0.500, -1.13, 1.400, -1.18, -1.59, 0.480, 0.730,
        0.950, 1.910, 0.520, 0.070, -1.27, -0.51, -0.21,
    ],
    [
        0.250, 0.040, -0.72, -0.62, 0.610, 0.160, -0.40, 0.730, -1.10, 0.530, 0.260, -0.64, -0.07,
        -0.69, -1.76, -0.26, -0.18, 0.540, 0.370, 0.020,
    ],
];
