// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/SEQUENCE/NeedlemanWunsch.h`: the global-alignment
//! similarity score [`ConsensusIDAlgorithmPEPMatrix`] uses.
//!
//! The header lives in `ANALYSIS/SEQUENCE`, which has no Rust domain module
//! yet; it sits here beside its only consumer.
//!
//! [`ConsensusIDAlgorithmPEPMatrix`]: crate::analysis::id::consensus_id::ConsensusIDAlgorithmPEPMatrix

use crate::{Error, Result};

/// Source `INT16_MAX`, the `J`/`O`/`U` rows and columns of both matrices.
const F: i32 = i16::MAX as i32;

/// Longest sequence [`NeedlemanWunsch::align`] accepts.
///
/// The source accumulates in a 32-bit `int`; with the largest matrix entry
/// (`INT16_MAX`) that cannot overflow below 65,536 residues. The ceiling also
/// bounds the rolling buffer, which is one row of the second sequence.
pub const MAX_ALIGNMENT_LENGTH: usize = 65_535;

/// Largest product of the two sequence lengths [`NeedlemanWunsch::align`]
/// evaluates: the dynamic-programming work in cells.
pub const MAX_ALIGNMENT_CELLS: usize = 100_000_000;

#[rustfmt::skip]
static MATRICES: [[[i32; 26]; 26]; 2] = [
    // identity
    [
        //  A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
        [1, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // A
        [0, 1, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // B
        [0, 0, 1, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // C
        [0, 0, 0, 1, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // D
        [0, 0, 0, 0, 1, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // E
        [0, 0, 0, 0, 0, 1, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // F
        [0, 0, 0, 0, 0, 0, 1, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // G
        [0, 0, 0, 0, 0, 0, 0, 1, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // H
        [0, 0, 0, 0, 0, 0, 0, 0, 1, F, 0, 1, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // I
        [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F], // J
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 1, 0, 0, 0, F, 0, 1, 0, 0, 0, F, 0, 0, 0, 0, 0], // K
        [0, 0, 0, 0, 0, 0, 0, 0, 1, F, 0, 1, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // L
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 1, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // M
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 1, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // N
        [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F], // O
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 1, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // P
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 1, 0, 0, 0, F, 0, 1, 0, 0, 0, F, 0, 0, 0, 0, 0], // Q
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 1, 0, 0, F, 0, 0, 0, 0, 0], // R
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 1, 0, F, 0, 0, 0, 0, 0], // S
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 1, F, 0, 0, 0, 0, 0], // T
        [F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F], // U
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 1, 0, 0, 0, 0], // V
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 1, 0, 0, 0], // W
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 0], // X
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 1, 0], // Y
        [0, 0, 0, 0, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, F, 0, 0, 0, 0, 0, F, 0, 0, 0, 0, 1], // Z
    ],
    // PAM30MS
    [
        //   A    B    C    D    E    F    G    H    I  J    K    L    M    N  O    P    Q    R    S    T  U    V    W  X    Y    Z
        [  6,  -7,  -6,  -3,  -2,  -8,  -2,  -7,  -5, F,  -7,  -6,  -5,  -4, F,  -2,  -4,  -7,   0,  -1, F,  -2, -13, 0,  -8,  -6], // A
        [ -7,   5, -11,  -7,  -7, -12,  -8,  -4,  -6, F,   5,  -7,  -3,  -4, F,  -5,  -3,   5,  -4,  -5, F,  -9,  -7, 0, -10,   1], // B
        [ -6, -11,  10, -14, -14, -13,  -9,  -7,  -6, F, -14, -11, -13, -11, F,  -8, -14,  -8,  -3,  -8, F,  -6, -15, 0,  -4, -14], // C
        [ -3,  -7, -14,   8,   2, -15,  -3,  -4,  -7, F,  -4, -10, -11,   2, F,  -8,  -2, -10,  -4,  -5, F,  -8, -15, 0, -11,  -3], // D
        [ -2,  -7, -14,   2,   8, -14,  -4,  -5,  -5, F,  -4,  -7,  -7,  -2, F,  -5,   1,  -9,  -4,  -6, F,  -6, -17, 0,  -8,  -2], // E
        [ -8, -12, -13, -15, -14,   9,  -9,  -6,  -2, F, -14,  -3,  -4,  -9, F, -10, -13,  -9,  -6,  -9, F,  -8,  -4, 0,   2, -14], // F
        [ -2,  -8,  -9,  -3,  -4,  -9,   6,  -9, -11, F,  -7, -11,  -8,  -3, F,  -6,  -7,  -9,  -2,  -6, F,  -5, -15, 0, -14,  -7], // G
        [ -7,  -4,  -7,  -4,  -5,  -6,  -9,   9,  -9, F,  -6,  -8, -10,   0, F,  -4,   1,  -2,  -6,  -7, F,  -6,  -7, 0,  -3,  -3], // H
        [ -5,  -6,  -6,  -7,  -5,  -2, -11,  -9,   8, F,  -6,   5,  -1,  -5, F,  -8,  -8,  -5,  -7,  -2, F,   2, -14, 0,  -6,  -7], // I
        [  F,   F,   F,   F,   F,   F,   F,   F,   F, F,   F,   F,   F,   F, F,   F,   F,   F,   F,   F, F,   F,   F, F,   F,   F], // J
        [ -7,   5, -14,  -4,  -4, -14,  -7,  -6,  -6, F,   7,  -7,  -2,  -1, F,  -6,  -3,   0,  -4,  -3, F,  -9, -12, 0,  -9,   4], // K
        [ -6,  -7, -11, -10,  -7,  -3, -11,  -8,   5, F,  -7,   5,   0,  -6, F,  -8,  -7,  -7,  -8,  -5, F,   0, -10, 0,  -7,  -7], // L
        [ -5,  -3, -13, -11,  -7,  -4,  -8, -10,  -1, F,  -2,   0,  11,  -9, F,  -8,  -4,  -4,  -5,  -4, F,  -1, -13, 0, -11,  -3], // M
        [ -4,  -4, -11,   2,  -2,  -9,  -3,   0,  -5, F,  -1,  -6,  -9,   8, F,  -6,  -3,  -6,   0,  -2, F,  -8,  -8, 0,  -4,  -2], // N
        [  F,   F,   F,   F,   F,   F,   F,   F,   F, F,   F,   F,   F,   F, F,   F,   F,   F,   F,   F, F,   F,   F, F,   F,   F], // O
        [ -2,  -5,  -8,  -8,  -5, -10,  -6,  -4,  -8, F,  -6,  -8,  -8,  -6, F,   8,  -3,  -4,  -2,  -4, F,  -6, -14, 0, -13,  -5], // P
        [ -4,  -3, -14,  -2,   1, -13,  -7,   1,  -8, F,  -3,  -7,  -4,  -3, F,  -3,   8,  -2,  -5,  -5, F,  -7, -13, 0, -12,   4], // Q
        [ -7,   5,  -8, -10,  -9,  -9,  -9,  -2,  -5, F,   0,  -7,  -4,  -6, F,  -4,  -2,   8,  -3,  -6, F,  -8,  -2, 0,  10,  -1], // R
        [  0,  -4,  -3,  -4,  -4,  -6,  -2,  -6,  -7, F,  -4,  -8,  -5,   0, F,  -2,  -5,  -3,   6,   0, F,  -6,  -5, 0,  -7,  -5], // S
        [ -1,  -5,  -8,  -5,  -6,  -9,  -6,  -7,  -2, F,  -3,  -5,  -4,  -2, F,  -4,  -5,  -6,   0,   7, F,  -3, -13, 0,  -6,  -4], // T
        [  F,   F,   F,   F,   F,   F,   F,   F,   F, F,   F,   F,   F,   F, F,   F,   F,   F,   F,   F, F,   F,   F, F,   F,   F], // U
        [ -2,  -9,  -6,  -8,  -6,  -8,  -5,  -6,   2, F,  -9,   0,  -1,  -8, F,  -6,  -7,  -8,  -6,  -3, F,   7, -15, 0,  -7,  -8], // V
        [-13,  -7, -15, -15, -17,  -4, -15,  -7, -14, F, -12, -10, -13,  -8, F, -14, -13,  -2,  -5, -13, F, -15,  13, 0,  -5, -13], // W
        [  0,   0,   0,   0,   0,   0,   0,   0,   0, F,   0,   0,   0,   0, F,   0,   0,   0,   0,   0, F,   0,   0, 0,   0,   0], // X
        [ -8, -10,  -4, -11,  -8,   2, -14,  -3,  -6, F,  -9,  -7, -11,  -4, F, -13, -12, -10,  -7,  -6, F,  -7,  -5, 0,  10, -11], // Y
        [ -6,   1, -14,  -3,  -2, -14,  -7,  -3,  -7, F,   4,  -7,  -3,  -2, F,  -5,   4,  -1,  -5,  -4, F,  -8, -13, 0, -11,   4], // Z
    ],
];

/// Substitution matrices bundled with [`NeedlemanWunsch`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum ScoringMatrix {
    /// +1 on the diagonal, 0 elsewhere, plus +1 for the isobaric pairs `I`/`L`
    /// and `K`/`Q`; `J`/`O`/`U` rows and columns hold `INT16_MAX`.
    Identity,
    /// PAM30 extended for MS-style substitutions (the default).
    #[default]
    PAM30MS,
}

impl ScoringMatrix {
    /// Source `NeedlemanWunsch::NamesOfScoringMatrices`, in enum order.
    pub const NAMES: [&'static str; 2] = ["identity", "PAM30MS"];

    /// The name [`ScoringMatrix::from_name`] accepts for this matrix.
    pub fn name(self) -> &'static str {
        Self::NAMES[self as usize]
    }

    /// The matrix with this exact (case-sensitive) name.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] for any name not in
    /// [`ScoringMatrix::NAMES`]; the message lists the valid names, as the
    /// source `Exception::IllegalArgument` does.
    pub fn from_name(name: &str) -> Result<Self> {
        match name {
            "identity" => Ok(Self::Identity),
            "PAM30MS" => Ok(Self::PAM30MS),
            _ => Err(Error::InvalidValue(format!(
                "Matrix is not known! Valid choices are: {}",
                Self::NAMES.join(", ")
            ))),
        }
    }

    fn table(self) -> &'static [[i32; 26]; 26] {
        &MATRICES[self as usize]
    }
}

/// Global-alignment similarity score for two amino-acid sequences using the
/// Needleman-Wunsch algorithm.
///
/// Computes the score (not the alignment itself) of the best global alignment
/// of two amino-acid sequences against a substitution matrix plus a linear gap
/// penalty, with a two-row rolling buffer.
///
/// # Notes
///
/// The source indexes the matrix by `c - 'A'`, so lowercase or non-ASCII input
/// reads out of bounds; [`NeedlemanWunsch::align`] refuses such input instead.
///
/// The source documentation says the `INT16_MAX` entries for `J`, `O` and `U`
/// "forbid" alignments involving them. They do the opposite: the recurrence
/// takes the maximum, so an `INT16_MAX` substitution is the most attractive
/// move available. The port reproduces the arithmetic, not the stated intent.
///
/// ```
/// use openms::analysis::id::needleman_wunsch::{NeedlemanWunsch, ScoringMatrix};
/// let aligner = NeedlemanWunsch::new(ScoringMatrix::Identity, 5);
/// assert_eq!(aligner.align("PEPTIDE", "PEPTIDE")?, 7);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct NeedlemanWunsch {
    /// The substitution matrix [`NeedlemanWunsch::align`] scores against.
    pub matrix: ScoringMatrix,
    /// Linear gap penalty, subtracted once per gap step. The source accepts any
    /// value; a negative one rewards gaps.
    pub penalty: i32,
}

impl Default for NeedlemanWunsch {
    /// [`ScoringMatrix::PAM30MS`] and a gap penalty of 5, as the source.
    fn default() -> Self {
        Self {
            matrix: ScoringMatrix::PAM30MS,
            penalty: 5,
        }
    }
}

impl NeedlemanWunsch {
    /// An aligner with this matrix and linear gap penalty.
    pub fn new(matrix: ScoringMatrix, penalty: i32) -> Self {
        Self { matrix, penalty }
    }

    /// Select the matrix by name (source `setMatrix(const std::string&)`).
    ///
    /// # Errors
    ///
    /// As [`ScoringMatrix::from_name`]; the aligner is unchanged on error.
    pub fn set_matrix_by_name(&mut self, name: &str) -> Result<()> {
        self.matrix = ScoringMatrix::from_name(name)?;
        Ok(())
    }

    /// Score of the best global alignment of `seq1` against `seq2`.
    ///
    /// The result is the source `int`, widened: both are exact for every input
    /// this accepts.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when either sequence contains anything
    /// but uppercase ASCII `A`-`Z` (source undefined behaviour), is longer than
    /// [`MAX_ALIGNMENT_LENGTH`], or when the product of the lengths exceeds
    /// [`MAX_ALIGNMENT_CELLS`].
    pub fn align(&self, seq1: &str, seq2: &str) -> Result<i64> {
        let a = residue_indices(seq1)?;
        let b = residue_indices(seq2)?;
        if a.len().saturating_mul(b.len()) > MAX_ALIGNMENT_CELLS {
            return Err(Error::InvalidValue(format!(
                "alignment of {} x {} residues exceeds {MAX_ALIGNMENT_CELLS} cells",
                a.len(),
                b.len()
            )));
        }
        let table = self.matrix.table();
        let gap = i64::from(self.penalty);
        let mut previous: Vec<i64> = (0..=b.len() as i64).map(|j| -j * gap).collect();
        let mut current = vec![0_i64; b.len() + 1];
        for (i, &x) in a.iter().enumerate() {
            current[0] = -((i as i64) + 1) * gap;
            for (j, &y) in b.iter().enumerate() {
                let substitution = previous[j] + i64::from(table[x][y]);
                current[j + 1] = (current[j] - gap)
                    .max(previous[j + 1] - gap)
                    .max(substitution);
            }
            std::mem::swap(&mut previous, &mut current);
        }
        Ok(previous[b.len()])
    }
}

fn residue_indices(sequence: &str) -> Result<Vec<usize>> {
    if sequence.len() > MAX_ALIGNMENT_LENGTH {
        return Err(Error::InvalidValue(format!(
            "alignment input of {} residues exceeds {MAX_ALIGNMENT_LENGTH}",
            sequence.len()
        )));
    }
    sequence
        .bytes()
        .map(|c| {
            if c.is_ascii_uppercase() {
                Ok(usize::from(c - b'A'))
            } else {
                Err(Error::InvalidValue(
                    "alignment input must contain only uppercase ASCII letters A-Z".into(),
                ))
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // Source `NeedlemanWunsch_test.cpp` at bc9cc12, lines 29-69.
    const SEQ1: &str = "IGGATLIGQLAIQQAHVHL";
    const SEQ2: &str = "IGGATLIGALDQVVAQQAHVHL";

    #[test]
    fn identity_matrix_class_test_values() {
        let aligner = NeedlemanWunsch::new(ScoringMatrix::Identity, 5);
        assert_eq!(aligner.align(SEQ1, SEQ2).unwrap(), 1);
        assert_eq!(aligner.align(SEQ1, SEQ1).unwrap(), 19);
        assert_eq!(aligner.align(SEQ2, SEQ2).unwrap(), 22);
    }

    #[test]
    fn pam30ms_matrix_class_test_values() {
        let mut aligner = NeedlemanWunsch::new(ScoringMatrix::Identity, 5);
        aligner.matrix = ScoringMatrix::PAM30MS;
        assert_eq!(aligner.align(SEQ1, SEQ2).unwrap(), 93);
        assert_eq!(aligner.align(SEQ1, SEQ1).unwrap(), 131);
        assert_eq!(aligner.align(SEQ2, SEQ2).unwrap(), 151);
    }

    #[test]
    fn matrix_by_name_class_test_values() {
        let mut aligner = NeedlemanWunsch::new(ScoringMatrix::PAM30MS, 5);
        assert!(aligner.set_matrix_by_name("Identity").is_err());
        assert_eq!(aligner.matrix, ScoringMatrix::PAM30MS);
        aligner.set_matrix_by_name("identity").unwrap();
        assert_eq!(aligner.align(SEQ1, SEQ2).unwrap(), 1);
        assert_eq!(aligner.align(SEQ1, SEQ1).unwrap(), 19);
        assert_eq!(aligner.align(SEQ2, SEQ2).unwrap(), 22);
    }

    #[test]
    fn penalty_class_test_value() {
        let mut aligner = NeedlemanWunsch::new(ScoringMatrix::PAM30MS, 5);
        aligner.penalty = 1;
        assert_eq!(aligner.align(SEQ1, SEQ2).unwrap(), 113);
        assert_eq!(aligner.penalty, 1);
    }

    #[test]
    fn default_is_pam30ms_with_penalty_five() {
        assert_eq!(
            NeedlemanWunsch::default(),
            NeedlemanWunsch::new(ScoringMatrix::PAM30MS, 5)
        );
    }

    #[test]
    fn refuses_input_the_source_reads_out_of_bounds() {
        let aligner = NeedlemanWunsch::default();
        assert!(aligner.align("PEPtIDE", "PEPTIDE").is_err());
        assert!(aligner.align("PEPTIDE", "PEP*").is_err());
    }

    #[test]
    fn empty_sequences_score_only_gaps() {
        // Derived from the recurrence: the first row is -j * penalty.
        let aligner = NeedlemanWunsch::default();
        assert_eq!(aligner.align("", "").unwrap(), 0);
        assert_eq!(aligner.align("", "PEP").unwrap(), -15);
        assert_eq!(aligner.align("PEP", "").unwrap(), -15);
    }
}
