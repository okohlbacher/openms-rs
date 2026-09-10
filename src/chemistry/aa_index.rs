// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Selected amino-acid indices and estimated peptide gas-phase basicity.
//!
//! Tables and ordinary arithmetic follow OpenMS4-core revision `7c029e8`,
//! [`AAIndex.h`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8/src/openms/include/OpenMS/CHEMISTRY/AAIndex.h).
//! These properties use parent residues and ignore every sequence annotation.

use super::AASequence;
use crate::{Error, Result};

/// The ten amino-acid scales exposed by the source AAIndex class.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AAIndexScale {
    /// Kerr-constant increments (Khanarian–Moore, 1980).
    Khag800101,
    /// Relative population of conformational state E (Vasquez et al., 1983).
    Vasm830103,
    /// Hydropathy, 36% accessibility (Naderi-Manesh et al., 2001).
    Nadh010106,
    /// Hydropathy, 50% accessibility (Naderi-Manesh et al., 2001).
    Nadh010107,
    /// RP-HPLC C8 hydrophobicity coefficient (Wilce et al., 1995).
    Wilm950102,
    /// Information measure for extended conformation without H-bond (Robson–Suzuki, 1976).
    Robb760107,
    /// Optimized average non-bonded energy per atom (Oobatake et al., 1985).
    Oobm850104,
    /// Positive-charge indicator (Fauchere et al., 1988).
    Fauj880111,
    /// Helix-coil equilibrium constant (Finkelstein–Ptitsyn, 1977).
    Fina770101,
    /// Signal-sequence helical potential (Argos et al., 1982).
    Argp820102,
}

impl AAIndexScale {
    /// Every supported scale in the source function order.
    pub const ALL: [Self; 10] = [
        Self::Khag800101,
        Self::Vasm830103,
        Self::Nadh010106,
        Self::Nadh010107,
        Self::Wilm950102,
        Self::Robb760107,
        Self::Oobm850104,
        Self::Fauj880111,
        Self::Fina770101,
        Self::Argp820102,
    ];

    /// The original AAindex accession, preserving its uppercase spelling.
    pub const fn accession(self) -> &'static str {
        match self {
            Self::Khag800101 => "KHAG800101",
            Self::Vasm830103 => "VASM830103",
            Self::Nadh010106 => "NADH010106",
            Self::Nadh010107 => "NADH010107",
            Self::Wilm950102 => "WILM950102",
            Self::Robb760107 => "ROBB760107",
            Self::Oobm850104 => "OOBM850104",
            Self::Fauj880111 => "FAUJ880111",
            Self::Fina770101 => "FINA770101",
            Self::Argp820102 => "ARGP820102",
        }
    }

    /// Look up one uppercase canonical residue. Unknown or lowercase codes error.
    pub fn value(self, residue: char) -> Result<f64> {
        Ok(INDEX_VALUES[self as usize][canonical_index(residue)?])
    }
}

/// Source residue indicators and peptide gas-phase basicity.
#[derive(Clone, Copy, Debug, Default)]
pub struct AAIndex;

impl AAIndex {
    /// Literal source indicator: A/G/F/I/M/L/P/V; every other character gives zero.
    pub fn aliphatic(residue: char) -> f64 {
        f64::from(matches!(
            residue,
            'A' | 'G' | 'F' | 'I' | 'M' | 'L' | 'P' | 'V'
        ))
    }

    /// Literal source indicator: D/E; every other character gives zero.
    pub fn acidic(residue: char) -> f64 {
        f64::from(matches!(residue, 'D' | 'E'))
    }

    /// Literal source indicator: K/R/H/W, including W; all other codes give zero.
    pub fn basic(residue: char) -> f64 {
        f64::from(matches!(residue, 'K' | 'R' | 'H' | 'W'))
    }

    /// Literal source indicator: S/T/Y/H/C/N/Q/W; all other codes give zero.
    pub fn polar(residue: char) -> f64 {
        f64::from(matches!(
            residue,
            'S' | 'T' | 'Y' | 'H' | 'C' | 'N' | 'Q' | 'W'
        ))
    }

    /// Estimated gas-phase basicity using the source's old gas constant and
    /// final division by ln(2). The source default temperature is 500 kelvin.
    ///
    /// Temperature must be finite and positive, and every parent residue must
    /// be canonical. All annotations are ignored. The first residue's sidechain
    /// is omitted, while later zero-valued sidechains still contribute exp(0).
    /// Ordinary finite arithmetic retains source split pairing and summation.
    /// Overflow retries an equivalent energy-domain log-sum-exp; an underflowed
    /// R*T uses its maximum-energy limit. Empty sequences use the exact one-site
    /// identity at all temperatures, avoiding the source's high-T cancellation.
    ///
    /// At most 1,000,000 residues and 50,000,000 work units are accepted. Worst-case
    /// work is charged before traversal; this calculation allocates no sequence
    /// or site-energy copies and never changes the input.
    pub fn calculate_gb(sequence: &AASequence, temperature_kelvin: f64) -> Result<f64> {
        if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
            return Err(invalid(
                "gas-basicity temperature must be finite and positive",
            ));
        }
        let n = sequence.len();
        // Cover validation and up to three complete site passes, including the
        // direct attempt before an overflow retry. Each split has <=3 lookups.
        let work = n.checked_add(1).and_then(|n| n.checked_mul(10));
        check_property_limits(n, work)?;
        let residues = sequence.as_str().as_bytes();
        for &residue in residues {
            canonical_index(char::from(residue))?;
        }
        let ln_two = 2.0_f64.ln();
        if residues.is_empty() {
            return finite((GB_N_TERM + GB_C_TERM) / ln_two);
        }
        let rt = GAS_CONSTANT_KJ * temperature_kelvin;
        if rt > 0.0 {
            let mut k_app = 0.0;
            for split in 0..=n {
                let (backbone, sidechain) = split_energies(residues, split)?;
                let mut contribution = (backbone / rt).exp();
                if let Some(sidechain) = sidechain {
                    contribution += (sidechain / rt).exp();
                }
                k_app += contribution;
                if !k_app.is_finite() {
                    break;
                }
            }
            let result = rt * k_app.ln() / ln_two;
            if result.is_finite() {
                return Ok(result);
            }
        }
        let mut maximum = 0.0_f64;
        for split in 0..=n {
            let (backbone, sidechain) = split_energies(residues, split)?;
            maximum = maximum.max(backbone);
            if let Some(sidechain) = sidechain {
                maximum = maximum.max(sidechain);
            }
        }
        if rt == 0.0 {
            // With bounded site energies/counts, the entropy correction is far
            // below a result ULP when this positive product rounds to zero.
            return finite(maximum / ln_two);
        }
        let mut scaled_sum = 0.0;
        for split in 0..=n {
            let (backbone, sidechain) = split_energies(residues, split)?;
            let mut contribution = ((backbone - maximum) / rt).exp();
            if let Some(sidechain) = sidechain {
                contribution += ((sidechain - maximum) / rt).exp();
            }
            scaled_sum += contribution;
        }
        finite((maximum + rt * scaled_sum.ln()) / ln_two)
    }
}

fn split_energies(residues: &[u8], split: usize) -> Result<(f64, Option<f64>)> {
    let left = if split == 0 {
        GB_N_TERM
    } else {
        GB_LEFT[canonical_index(char::from(residues[split - 1]))?]
    };
    if split == residues.len() {
        return Ok((left + GB_C_TERM, None));
    }
    let right = canonical_index(char::from(residues[split]))?;
    Ok((
        left + GB_RIGHT[right],
        (split > 0).then_some(GB_SIDECHAIN[right]),
    ))
}

pub(super) const MAX_PROPERTY_RESIDUES: usize = 1_000_000;
pub(super) const MAX_PROPERTY_WORK: usize = 50_000_000;

pub(super) fn check_property_limits(residues: usize, work: Option<usize>) -> Result<()> {
    if residues > MAX_PROPERTY_RESIDUES {
        return Err(invalid("peptide property input exceeds 1,000,000 residues"));
    }
    if work.is_none_or(|work| work > MAX_PROPERTY_WORK) {
        return Err(invalid(
            "peptide property calculation exceeds 50,000,000 work units",
        ));
    }
    Ok(())
}

pub(super) fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

pub(super) fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("peptide property result must be finite"))
    }
}

// Shared canonical column order: ACDEFGHIKLMNPQRSTVWY. The permissive source
// indicator functions intentionally do not use this checked numeric lookup.
pub(super) fn canonical_index(residue: char) -> Result<usize> {
    match residue {
        'A' => Ok(0),
        'C' => Ok(1),
        'D' => Ok(2),
        'E' => Ok(3),
        'F' => Ok(4),
        'G' => Ok(5),
        'H' => Ok(6),
        'I' => Ok(7),
        'K' => Ok(8),
        'L' => Ok(9),
        'M' => Ok(10),
        'N' => Ok(11),
        'P' => Ok(12),
        'Q' => Ok(13),
        'R' => Ok(14),
        'S' => Ok(15),
        'T' => Ok(16),
        'V' => Ok(17),
        'W' => Ok(18),
        'Y' => Ok(19),
        _ => Err(Error::InvalidValue(format!(
            "no peptide property value for residue {residue:?}; expected an uppercase canonical residue"
        ))),
    }
}

// Exact pinned AAIndex.h values. Columns use the canonical order above.
const INDEX_VALUES: [[f64; 20]; 10] = [
    // KHAG800101
    [
        49.1, 0.0, 0.0, 0.0, 54.7, 64.6, 75.7, 18.9, 0.0, 15.6, 6.8, -3.6, 43.8, 20.0, 133.0, 44.4,
        31.0, 29.5, 70.5, 0.0,
    ],
    // VASM830103
    [
        0.159, 0.187, 0.283, 0.206, 0.682, 0.049, 0.233, 0.581, 0.159, 0.083, 0.198, 0.385, 0.366,
        0.236, 0.194, 0.150, 0.074, 0.301, 0.463, 0.737,
    ],
    // NADH010106
    [
        5.0, 224.0, 45.0, -8.0, 117.0, -47.0, -50.0, 83.0, -38.0, 82.0, 83.0, -77.0, -103.0, -67.0,
        -57.0, -41.0, 79.0, 117.0, 130.0, 27.0,
    ],
    // NADH010107
    [
        -2.0, 329.0, 248.0, 117.0, 120.0, -66.0, -70.0, 28.0, 115.0, 36.0, 62.0, -97.0, -132.0,
        -37.0, -41.0, -52.0, 174.0, 114.0, 179.0, -7.0,
    ],
    // WILM950102
    [
        2.62, 0.73, -2.84, -0.45, 9.14, -1.15, -0.74, 4.38, -2.78, 6.57, -3.12, -1.27, -0.12,
        -1.69, 1.26, -1.39, 1.81, 2.30, 5.91, 1.39,
    ],
    // ROBB760107
    [
        0.0, 5.4, -2.6, 3.1, 0.7, -3.4, 0.8, -0.1, -3.1, -3.7, -2.1, -2.0, 7.4, 2.4, 1.1, 1.3, 0.0,
        2.7, -3.4, 4.8,
    ],
    // OOBM850104
    [
        -2.49, -3.13, 8.86, 4.04, -6.64, -0.56, 4.22, -10.87, -9.97, -7.16, -4.96, 2.27, 5.19,
        1.79, 2.55, -1.60, -4.75, -3.97, -17.84, 9.25,
    ],
    // FAUJ880111
    [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0,
    ],
    // FINA770101
    [
        1.08, 0.95, 0.85, 1.15, 1.10, 0.55, 1.00, 1.05, 1.15, 1.25, 1.15, 0.85, 0.71, 0.95, 1.05,
        0.75, 0.75, 0.95, 1.10, 1.10,
    ],
    // ARGP820102
    [
        1.18, 1.89, 0.05, 0.11, 1.96, 0.49, 0.31, 1.45, 0.06, 3.23, 2.67, 0.23, 0.76, 0.72, 0.20,
        0.97, 0.84, 1.08, 0.77, 0.39,
    ],
];

const GB_LEFT: [f64; 20] = [
    881.82, 881.15, 880.02, 880.10, 881.08, 881.17, 881.27, 880.99, 880.06, 881.88, 881.38, 881.18,
    881.25, 881.50, 882.98, 881.08, 881.14, 881.17, 881.31, 881.20,
];
// The tabulated arginine energy delta 6.28 is not the mathematical constant tau.
#[allow(clippy::approx_constant)]
const GB_RIGHT: [f64; 20] = [
    0.0, -0.69, -0.63, -0.39, 0.03, 0.92, -0.19, -1.17, -0.71, -0.09, 0.30, 1.56, 11.75, 4.10,
    6.28, 0.98, 1.21, -0.90, 0.10, -0.38,
];
const GB_SIDECHAIN: [f64; 20] = [
    0.0, 0.0, 784.0, 790.0, 0.0, 0.0, 927.84, 0.0, 926.74, 0.0, 830.0, 864.94, 0.0, 865.25, 1000.0,
    775.0, 780.0, 0.0, 909.53, 790.0,
];
const GB_N_TERM: f64 = 916.84;
const GB_C_TERM: f64 = -95.82;
// Deliberately retain the old Constants.h multiplication/division order.
const GAS_CONSTANT_KJ: f64 = (6.022_136_7e23 * 1.380_657e-23) / 1000.0;
