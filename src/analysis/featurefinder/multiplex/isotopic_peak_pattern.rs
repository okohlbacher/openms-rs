// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The pattern of isotopic peaks a labelled peptide multiplet forms.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexIsotopicPeakPattern.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use super::delta_masses::MultiplexDeltaMasses;
use crate::concept::constants::C13C12_MASSDIFF_U;
use crate::{Error, Result};

/// Ceiling on the number of m/z shifts one pattern may hold
/// (`peaks_per_peptide` times the number of mass shifts), checked before the
/// shift list is allocated. The source has no such bound.
pub const MAX_MZ_SHIFTS: usize = 1_000_000;

/// Pattern of isotopic peaks of a peptide multiplet (source
/// `MultiplexIsotopicPeakPattern`).
///
/// Groups of peptides appear as characteristic patterns of isotopic peaks in
/// MS1 spectra. For an Arg6-labelled SILAC peptide pair of charge 2+ with three
/// isotopic peaks we expect peaks at relative m/z shifts of 0, 0.5, 1, 3, 3.5
/// and 4 Th.
#[derive(Clone, Debug, PartialEq)]
pub struct MultiplexIsotopicPeakPattern {
    /// m/z shifts between isotopic peaks, `peaks_per_peptide` for each mass
    /// shift, peptide by peptide (source `mz_shifts_`).
    mz_shifts: Vec<f64>,
    /// Charge (source `charge_`).
    charge: i32,
    /// Number of isotopic peaks in each peptide (source `peaks_per_peptide_`).
    peaks_per_peptide: i32,
    /// Mass shifts between peptides, including the zero mass shift of the first
    /// peptide (source `mass_shifts_`).
    mass_shifts: MultiplexDeltaMasses,
    /// Index in the mass shift list (source `mass_shift_index_`).
    mass_shift_index: i32,
}

impl MultiplexIsotopicPeakPattern {
    /// A pattern of charge `charge` with `peaks_per_peptide` isotopic peaks per
    /// peptide, the mass shifts `mass_shifts` and the index `mass_shift_index`
    /// in the mass shift list (source constructor).
    ///
    /// The m/z shift of isotope `j` of peptide `i` is `(delta_mass_i + j *
    /// C13C12_MASSDIFF_U) / charge`, evaluated in that order. A negative
    /// `peaks_per_peptide` yields no m/z shifts, as the source's loop does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for a zero charge, where the source
    /// divides by zero and stores infinities and NaNs, and
    /// [`Error::InvalidRange`] when the pattern would hold more than
    /// [`MAX_MZ_SHIFTS`] m/z shifts.
    pub fn new(
        charge: i32,
        peaks_per_peptide: i32,
        mass_shifts: MultiplexDeltaMasses,
        mass_shift_index: i32,
    ) -> Result<Self> {
        if charge == 0 {
            return Err(Error::InvalidValue(
                "MultiplexIsotopicPeakPattern: charge must not be zero".into(),
            ));
        }
        let per_peptide = usize::try_from(peaks_per_peptide.max(0)).unwrap_or(0);
        let count = per_peptide
            .checked_mul(mass_shifts.delta_masses.len())
            .filter(|&n| n <= MAX_MZ_SHIFTS)
            .ok_or_else(|| {
                Error::InvalidRange(format!(
                    "MultiplexIsotopicPeakPattern: more than {MAX_MZ_SHIFTS} m/z shifts"
                ))
            })?;
        let mut mz_shifts = Vec::with_capacity(count);
        for delta in &mass_shifts.delta_masses {
            for j in 0..per_peptide {
                mz_shifts
                    .push((delta.delta_mass + j as f64 * C13C12_MASSDIFF_U) / f64::from(charge));
            }
        }
        Ok(Self {
            mz_shifts,
            charge,
            peaks_per_peptide,
            mass_shifts,
            mass_shift_index,
        })
    }

    /// The charge (source `getCharge`).
    pub fn charge(&self) -> i32 {
        self.charge
    }

    /// The number of isotopic peaks per peptide (source `getPeaksPerPeptide`).
    pub fn peaks_per_peptide(&self) -> i32 {
        self.peaks_per_peptide
    }

    /// The mass shifts (source `getMassShifts`, which returns a copy).
    pub fn mass_shifts(&self) -> &MultiplexDeltaMasses {
        &self.mass_shifts
    }

    /// The index in the mass shift list (source `getMassShiftIndex`).
    pub fn mass_shift_index(&self) -> i32 {
        self.mass_shift_index
    }

    /// The number of mass shifts, i.e. the number of peptides in the multiplet
    /// (source `getMassShiftCount`).
    pub fn mass_shift_count(&self) -> usize {
        self.mass_shifts.delta_masses.len()
    }

    /// The mass shift of peptide `i`, or `None` past the end (source
    /// `getMassShiftAt`, undefined behaviour past the end).
    pub fn mass_shift_at(&self, i: usize) -> Option<f64> {
        self.mass_shifts.delta_masses.get(i).map(|d| d.delta_mass)
    }

    /// The m/z shift at position `i`, or `None` past the end (source
    /// `getMZShiftAt`, undefined behaviour past the end). Position
    /// `peptide * peaks_per_peptide + isotope` addresses one isotope of one
    /// peptide.
    pub fn mz_shift_at(&self, i: usize) -> Option<f64> {
        self.mz_shifts.get(i).copied()
    }

    /// The number of m/z shifts (source `getMZShiftCount`).
    pub fn mz_shift_count(&self) -> usize {
        self.mz_shifts.len()
    }

    /// All m/z shifts in position order.
    pub fn mz_shifts(&self) -> &[f64] {
        &self.mz_shifts
    }
}
