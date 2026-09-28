// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Port of `ANALYSIS/ID/CometModification.h`: one Comet `variable_modXX`
//! entry and the rules for collapsing compatible entries into one line.
//!
//! Comet accepts several residues in one modification entry (for example `KR`
//! for lysine and arginine), which makes searches markedly faster. The merging
//! rules, as the source's brief states and its code implements:
//!
//! * two entries need the same mass (within
//!   [`CometModification::MASS_TOLERANCE`]) and the same binary group;
//! * plain amino-acid entries (`term_distance == -1`) merge with each other
//!   (`K` + `R` gives `KR`);
//! * a peptide N-term entry (`nc_term == 2`, `term_distance == 0`) merges with
//!   amino-acid entries (`n` + `RST` gives `nRST`), and likewise a peptide
//!   C-term entry (`nc_term == 3`; `c` + `KR` gives `cKR`);
//! * a protein-terminal entry (`nc_term` 0 or 1 with `term_distance == 0`)
//!   merges only with another protein-terminal entry of the same kind;
//! * peptide N-term and peptide C-term entries never merge with each other.
//!
//! When a terminal entry and an amino-acid entry merge, the terminal
//! specificity moves into the residue string as `n` or `c`, and
//! `term_distance` / `nc_term` become `-1` / `0`, following Comet's own
//! example `42.010565 nK 0 3 -1 0 0 0.0`. The mapping and the one known source
//! defect are recorded in `docs/COMET_MODIFICATION_SUPPORT.md`.

use crate::chemistry::{ResidueModification, TermSpecificity};
use crate::concept::fuzzy_string_comparator::format_g;

/// One Comet `variable_modXX` entry (source struct `CometModification`).
///
/// The fields mirror the `variable_modXX` columns. `Default` is the source's
/// default constructor: a zeroed entry with mass `0`, no residues, binary group
/// `0`, at most 5 per peptide, no terminal constraint, not required.
#[derive(Clone, Debug, PartialEq)]
pub struct CometModification {
    /// Modification mass difference, in Da.
    pub mass: f64,
    /// Residue(s) the modification applies to, e.g. `K`, `KR`, `n`, `nKR`.
    pub residues: String,
    /// Comet binary modification group (used e.g. for SILAC).
    pub binary_group: i32,
    /// Maximum number of occurrences of this modification per peptide.
    pub max_mods_per_peptide: i32,
    /// Terminal distance constraint: `-1` for none, `0` for terminal only.
    pub term_distance: i32,
    /// Terminal specificity: `0` protein N-term, `1` protein C-term, `2`
    /// peptide N-term, `3` peptide C-term. Only meaningful when
    /// `term_distance == 0`.
    pub nc_term: i32,
    /// Whether the modification must appear in every peptide.
    pub required: bool,
}

impl Default for CometModification {
    fn default() -> Self {
        Self {
            mass: 0.0,
            residues: String::new(),
            binary_group: 0,
            max_mods_per_peptide: 5,
            term_distance: -1,
            nc_term: 0,
            required: false,
        }
    }
}

impl CometModification {
    /// Absolute tolerance, in Da, that [`CometModification::is_mergeable_with`]
    /// allows between two masses.
    pub const MASS_TOLERANCE: f64 = 1e-6;

    /// Build an entry from an OpenMS modification (source constructor from a
    /// `const ResidueModification*`).
    ///
    /// Copies the modification's monoisotopic mass difference into `mass` and
    /// its origin residue into `residues`, then sets `term_distance` and
    /// `nc_term` from the term specificity:
    ///
    /// * peptide C-term: `0` and `3`;
    /// * peptide N-term: `0` and `2`;
    /// * protein N-term: `0` and `0`;
    /// * protein C-term: `0` and `1`;
    /// * anywhere: left at `-1` and `0`.
    ///
    /// A terminal modification whose origin is `X` (any residue) gets the
    /// residue string `n` or `c`, which Comet recognises; a terminal
    /// modification on a specific residue keeps that residue. This registry
    /// stores the origin of an any-residue terminal modification as `None`
    /// where the source stores `X`; both are treated as `X`, and an anywhere
    /// modification without an origin becomes `X`, as the source's character
    /// would be. `required` stays `false`.
    ///
    /// # Arguments
    ///
    /// * `modification` - the source modification; the source takes a pointer
    ///   and does not check it for null, a reference cannot be null
    /// * `binary_grp` - value for `binary_group`
    /// * `max_mods` - value for `max_mods_per_peptide`
    pub fn from_modification(
        modification: &ResidueModification,
        binary_grp: i32,
        max_mods: i32,
    ) -> Self {
        let origin = modification.origin().unwrap_or('X');
        let any_residue = origin == 'X';
        let mut entry = Self {
            mass: modification.diff_mono_mass(),
            residues: origin.to_string(),
            binary_group: binary_grp,
            max_mods_per_peptide: max_mods,
            ..Self::default()
        };
        let (terminal, nc_term) = match modification.term_specificity() {
            TermSpecificity::CTerm => ("c", 3),
            TermSpecificity::NTerm => ("n", 2),
            TermSpecificity::ProteinNTerm => ("n", 0),
            TermSpecificity::ProteinCTerm => ("c", 1),
            TermSpecificity::Anywhere => return entry,
        };
        if any_residue {
            entry.residues = terminal.into();
        }
        entry.term_distance = 0;
        entry.nc_term = nc_term;
        entry
    }

    fn is_protein_terminal(&self) -> bool {
        self.term_distance == 0 && (self.nc_term == 0 || self.nc_term == 1)
    }

    fn is_peptide_terminal(&self, nc_term: i32) -> bool {
        self.term_distance == 0 && self.nc_term == nc_term
    }

    /// Whether two entries can be combined into one Comet `variable_modXX`
    /// line, by the rules in the module documentation.
    ///
    /// `false` when the masses differ by more than
    /// [`CometModification::MASS_TOLERANCE`], when the binary groups differ,
    /// when only one of the two is a protein-terminal entry or they are
    /// different protein-terminal kinds, or when one is a peptide N-term and
    /// the other a peptide C-term entry. The mass test is the source's
    /// `|a - b| > tolerance` refusal, so a NaN mass, which fails that test,
    /// is not refused by it.
    pub fn is_mergeable_with(&self, other: &CometModification) -> bool {
        if (self.mass - other.mass).abs() > Self::MASS_TOLERANCE {
            return false;
        }
        if self.binary_group != other.binary_group {
            return false;
        }
        if self.is_protein_terminal() || other.is_protein_terminal() {
            return self.nc_term == other.nc_term && self.term_distance == other.term_distance;
        }
        let crossed = (self.is_peptide_terminal(2) && other.is_peptide_terminal(3))
            || (self.is_peptide_terminal(3) && other.is_peptide_terminal(2));
        !crossed
    }

    /// Merge a compatible entry into this one.
    ///
    /// Appends `other`'s residue characters that are not yet present, keeps
    /// the larger `max_mods_per_peptide` and ORs `required`. When exactly one
    /// of the two is terminal (`term_distance == 0`), the terminal information
    /// is already in the residue string, so `term_distance` / `nc_term` become
    /// `-1` / `0`; two terminal or two non-terminal entries keep both fields.
    ///
    /// Behaviour is only meaningful for pairs that
    /// [`CometModification::is_mergeable_with`] accepts; as in the source, the
    /// caller is responsible for that check, and this method does not repeat
    /// it.
    pub fn merge(&mut self, other: &CometModification) {
        for residue in other.residues.chars() {
            if !self.residues.contains(residue) {
                self.residues.push(residue);
            }
        }
        self.max_mods_per_peptide = self.max_mods_per_peptide.max(other.max_mods_per_peptide);
        let this_is_terminal = self.term_distance == 0;
        let other_is_terminal = other.term_distance == 0;
        if this_is_terminal != other_is_terminal {
            self.term_distance = -1;
            self.nc_term = 0;
        }
        self.required = self.required || other.required;
    }

    /// Render this entry as a Comet `variable_modXX` parameter line.
    ///
    /// ```text
    /// variable_mod<II> = <mass> <residues> <binary_group> <max_mods_per_peptide> <term_distance> <nc_term> <required> 0.0
    /// ```
    ///
    /// `II` is `index` zero-padded to two digits (wider indices print in full,
    /// as `std::setw` does). The mass uses the `std::ostream` default of six
    /// significant digits; `required` prints as `1` or `0`. The last column,
    /// the neutral loss, is always `0.0`, as OpenMS does not surface it.
    ///
    /// ```
    /// use openms::analysis::id::comet_modification::CometModification;
    ///
    /// let entry = CometModification {
    ///     mass: 42.010565,
    ///     residues: "nK".into(),
    ///     max_mods_per_peptide: 3,
    ///     ..CometModification::default()
    /// };
    /// assert_eq!(entry.to_comet_string(1), "variable_mod01 = 42.0106 nK 0 3 -1 0 0 0.0");
    /// ```
    pub fn to_comet_string(&self, index: usize) -> String {
        format!(
            "variable_mod{index:02} = {} {} {} {} {} {} {} 0.0",
            format_g(self.mass),
            self.residues,
            self.binary_group,
            self.max_mods_per_peptide,
            self.term_distance,
            self.nc_term,
            u8::from(self.required),
        )
    }

    /// Greedy first-fit merge of a list of entries (source
    /// `mergeModifications`).
    ///
    /// Walks `mods` in order; each entry is merged into the first compatible
    /// entry already in the result, or appended unchanged when there is none.
    /// The kept entries follow the input order and nothing is dropped. The
    /// cost is quadratic in the number of entries in the worst case; Comet
    /// modification lists are short.
    pub fn merge_modifications(mods: &[CometModification]) -> Vec<CometModification> {
        let mut merged: Vec<CometModification> = Vec::with_capacity(mods.len());
        for modification in mods {
            match merged
                .iter_mut()
                .find(|existing| existing.is_mergeable_with(modification))
            {
                Some(existing) => existing.merge(modification),
                None => merged.push(modification.clone()),
            }
        }
        merged
    }
}
