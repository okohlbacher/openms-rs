// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS contributors $
//! Molecular formulas, modified peptides, isotope patterns, and digestion.
//!
//! Mass tables and peptide chemistry follow OpenMS4-core revision `7c029e8`:
//! [`ElementDB.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8/src/openms/source/CHEMISTRY/ElementDB.cpp),
//! `EmpiricalFormula.cpp`, `ResidueDB.cpp`, `Residue.h`, and `ProteaseDigestion.cpp`.
//! All 84 declared element tables are included. Iridium uses the declared iridium
//! data, correcting an upstream call that accidentally initializes it as rhenium.
//!
//! Formula charge means addition/removal of **protons**, as in OpenMS; it is not
//! an electron-only ionization adjustment. The embedded modification registry
//! includes the pinned UniMod and OpenMS custom tables. B/Z/X preserve unresolved
//! chemistry, and numeric mass tags retain known mass without inventing a formula.
//! Fine isotope structure uses bounded native configuration enumeration. Arbitrary
//! enzyme regular expressions remain outside the implemented surface.

use crate::error::{Error, Result};
use std::collections::BTreeMap;
use std::fmt;
use std::str::FromStr;

pub mod aa_index;
pub mod adduct_info;
pub mod decoy_generator;
pub mod digestion;
mod elements;
pub mod hydrophobicity;
pub mod ims_alphabet;
pub mod ims_element;
pub mod ims_isotope_distribution;
pub use ims_alphabet::{IMSAlphabet, IMSAlphabetParser, IMSAlphabetTextParser};
pub use ims_element::IMSElement;
pub use ims_isotope_distribution::{IMSIsotopeDistribution, IMSIsotopeOptions, IMSIsotopePeak};
pub mod ims_mass_decomposer;
pub use ims_mass_decomposer::{IMSIntegerMassDecomposer, IMSMassDecomposer, IMSRealMassDecomposer};
pub mod ims_weights;
pub mod ion_naming;
pub use ims_weights::IMSWeights;
pub mod isoelectric_point;
pub mod mass_decomposition;
pub mod mass_decomposition_algorithm;
pub use mass_decomposition::MassDecomposition;
pub use mass_decomposition_algorithm::{
    DecompositionResidueSet, MassDecompositionAlgorithm, MassDecompositionOptions,
};
pub mod modified_na_sequence_generator;
pub mod na_sequence;
pub mod nucleic_acid_spectrum_generator;
pub mod ribonucleotide;
pub mod ribonucleotide_db;
pub mod rnase;
pub mod sequence_coverage;
pub use sequence_coverage::SequenceCoverage;
pub mod spectrum_annotator;
pub mod tagger;
pub use aa_index::{AAIndex, AAIndexScale};
pub use adduct_info::AdductInfo;
pub use decoy_generator::DecoyGenerator;
pub use digestion::{
    DigestedPeptide, DigestionEnzymeProtein, DigestionProduct, DigestionSpecificity,
    ProductValidation, Protease, ProteaseDB, ProteaseDigestion,
};
pub use hydrophobicity::{HydrophobicityProfile, HydrophobicityScale};
pub use isoelectric_point::{IsoelectricPoint, ProteomicsPkaScale};
pub use modified_na_sequence_generator::ModifiedNASequenceGenerator;
pub use na_sequence::{NAFragmentType, NASequence};
pub use nucleic_acid_spectrum_generator::NucleicAcidSpectrumGenerator;
pub use ribonucleotide::{Ribonucleotide, RibonucleotideRecord, RibonucleotideTermSpecificity};
pub use ribonucleotide_db::{
    RibonucleotideDB, RibonucleotideDiagnostic, RibonucleotideEntry, RibonucleotideLoadReport,
};
pub use rnase::{
    DigestedOligo, DigestionEnzymeRNA, DigestionEnzymeRNARecord, RNaseDB, RNaseDigestion,
};
pub use spectrum_annotator::SpectrumAnnotator;
pub use tagger::{Tagger, TaggerOptions};
pub mod fine_isotopes;
pub mod isotopes;
pub use fine_isotopes::{
    FineIsotopeConfiguration, FineIsotopeIterator, FineIsotopePatternGenerator, FineIsotopeStop,
};
pub use isotopes::{
    AveragineComposition, CoarseIsotopePatternGenerator, CoarseMassMode, FormulaEstimate,
    IsotopeDistribution, IsotopePeak,
};
pub mod cross_links;
pub mod modification_definitions;
mod modifications;
pub mod modified_peptides;
pub mod protein_cross_link;
mod sequence;
pub mod theoretical;
pub mod theoretical_xlms;
pub use cross_links::CrossLinksDB;
pub use modification_definitions::{
    ModificationDefinition, ModificationDefinitionsSet, ModificationMassMode, ModificationMatch,
    ModificationMatchOptions,
};
pub use modifications::{
    ModificationProvenance, ModificationRecord, ModificationsDB, NeutralLoss, OboLoadReport,
    OboReadOptions, ResidueModification, TermSpecificity,
};
pub use modified_peptides::ModifiedPeptideGenerator;
pub use protein_cross_link::{ProteinProteinCrossLink, ProteinProteinCrossLinkType};
pub use theoretical_xlms::{LossIndex, TheoreticalSpectrumGeneratorXLMS, XLMSLimits, XLMSOptions};
pub mod monosaccharide_db;
pub use monosaccharide_db::{Monosaccharide, MonosaccharideDB};
pub use sequence::{AASequence, MassTag, PeptideFragmentType, SequenceModification};
pub use theoretical::{
    TheoreticalIonIntensities, TheoreticalIonSeries, TheoreticalIsotopeModel,
    TheoreticalSpectrumGenerator,
};

/// Existing chemistry names retain the source constants and exact values.
pub use crate::constants::{C13C12_MASSDIFF_U, ELECTRON_MASS_U, PROTON_MASS_U};

/// An isotope's exact upstream tabulated mass and fractional natural abundance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Isotope {
    pub mass_number: u16,
    pub mass: f64,
    pub abundance: f64,
}

/// An immutable element with naturally occurring isotope data.
///
/// Obtain validated elements through [`element`] or [`element_table`]. Their
/// fields cannot be replaced, even on an owned copy:
///
/// ```compile_fail
/// use openms::chemistry::element;
/// let mut carbon = *element("C").unwrap();
/// carbon.isotopes = &[];
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Element {
    name: &'static str,
    symbol: &'static str,
    atomic_number: u8,
    isotopes: &'static [Isotope],
}

impl Element {
    pub fn name(&self) -> &'static str {
        self.name
    }

    pub fn symbol(&self) -> &'static str {
        self.symbol
    }

    pub fn atomic_number(&self) -> u8 {
        self.atomic_number
    }

    /// The nonempty isotope table for this element, in mass-number order.
    pub fn isotopes(&self) -> &'static [Isotope] {
        self.isotopes
    }

    /// Mass of the most abundant isotope, matching OpenMS's definition.
    pub fn mono_mass(&self) -> f64 {
        let mut most_abundant = &self.isotopes[0];
        for isotope in &self.isotopes[1..] {
            if isotope.abundance > most_abundant.abundance {
                most_abundant = isotope;
            }
        }
        most_abundant.mass
    }

    /// Abundance-weighted average mass; the tabulated abundances are not rescaled.
    pub fn average_mass(&self) -> f64 {
        self.isotopes.iter().map(|i| i.mass * i.abundance).sum()
    }
}

/// All elements declared by the pinned upstream ElementDB, in atomic-number order.
pub fn element_table() -> &'static [Element] {
    elements::ELEMENTS
}

/// Look up a natural element by its case-sensitive symbol or full English name.
pub fn element(symbol_or_name: &str) -> Option<&'static Element> {
    element_table()
        .iter()
        .find(|e| e.symbol == symbol_or_name || e.name == symbol_or_name)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Atom {
    symbol: &'static str,
    isotope: Option<u16>,
}

impl Atom {
    fn resolve(symbol: &str, isotope: Option<u16>) -> Result<Self> {
        let (symbol, isotope) = match (symbol, isotope) {
            ("D", None) => ("H", Some(2)),
            ("T", None) => ("H", Some(3)),
            _ => (symbol, isotope),
        };
        let entry = element(symbol)
            .filter(|e| e.symbol == symbol)
            .ok_or_else(|| Error::InvalidValue(format!("unknown element symbol {symbol:?}")))?;
        if let Some(number) = isotope {
            if !entry.isotopes.iter().any(|i| i.mass_number == number) {
                return Err(Error::InvalidValue(format!(
                    "isotope ({number}){symbol} is not in the element table"
                )));
            }
        }
        Ok(Self {
            symbol: entry.symbol,
            isotope,
        })
    }

    fn mass(self, average: bool) -> f64 {
        // Atoms can only be constructed from validated entries in the static table.
        let entry = element(self.symbol).expect("validated element");
        match self.isotope {
            Some(number) => {
                entry
                    .isotopes
                    .iter()
                    .find(|i| i.mass_number == number)
                    .expect("validated isotope")
                    .mass
            }
            None if average => entry.average_mass(),
            None => entry.mono_mass(),
        }
    }
}

/// A molecular composition with signed atom counts and a protonation charge.
///
/// Supported notation: `C6H12O6`, `(13)C6H12O6`, `D2O`, `C-1H2`, `C6H12O6+2`,
/// `C6H12O6-2`, and charge-only `+2`. A minus immediately after a symbol starts a
/// negative atom count: `H-2` removes two hydrogen atoms; `H1-2` has one hydrogen
/// atom and charge -2. Use [`Self::with_charge`] to avoid notation ambiguity.
/// Outer whitespace is trimmed. Group multipliers and interior whitespace are
/// rejected. Display writes explicit counts so negative charges round-trip.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EmpiricalFormula {
    atoms: BTreeMap<Atom, i32>,
    charge: i32,
}

impl EmpiricalFormula {
    pub fn parse(input: &str) -> Result<Self> {
        input.parse()
    }

    pub fn charge(&self) -> i32 {
        self.charge
    }

    pub fn with_charge(mut self, charge: i32) -> Self {
        self.charge = charge;
        self
    }

    /// True when there are no atoms, independently of charge.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Signed total atom count; isotope labels are counted separately.
    pub fn atom_count(&self) -> i64 {
        self.atoms.values().map(|&n| i64::from(n)).sum()
    }

    /// Number of stored keys for bounded graph measurement; no formula formatting.
    pub(crate) fn stored_atom_types(&self) -> usize {
        self.atoms.len()
    }

    /// Deterministic structural ordering for identification graph keys.
    pub(crate) fn graph_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.atoms
            .len()
            .cmp(&other.atoms.len())
            .then_with(|| self.charge.cmp(&other.charge))
            .then_with(|| self.atoms.iter().cmp(other.atoms.iter()))
    }

    /// Count an element or labeled isotope (`C`, `(13)C`, `D`); absent means zero.
    /// Natural `C` does not include explicitly labeled `(13)C` atoms.
    pub fn count(&self, symbol: &str) -> Result<i32> {
        let mut position = 0;
        let atom = parse_atom(symbol.as_bytes(), &mut position)?;
        if position != symbol.len() {
            return Err(Error::InvalidValue("expected one element symbol".into()));
        }
        Ok(self.atoms.get(&atom).copied().unwrap_or(0))
    }

    /// Sum of most-abundant isotope masses and `charge * PROTON_MASS_U`.
    pub fn mono_mass(&self) -> f64 {
        self.mass(false)
    }

    /// Sum of abundance-weighted masses and `charge * PROTON_MASS_U`.
    pub fn average_mass(&self) -> f64 {
        self.mass(true)
    }

    fn mass(&self, average: bool) -> f64 {
        self.atoms.iter().fold(
            f64::from(self.charge) * PROTON_MASS_U,
            |mass, (&atom, &count)| mass + atom.mass(average) * f64::from(count),
        )
    }

    /// Nonnegative mass-to-charge ratio for positive or negative protonation.
    /// Returns an error for zero charge or a negative calculated ion mass.
    pub fn mz(&self) -> Result<f64> {
        ion_mz(self.mono_mass(), self.charge)
    }

    /// Check that each atom count in `other` is no greater than this count.
    /// Missing atoms count as zero; signed counts are compared and charge is ignored.
    pub fn contains(&self, other: &Self) -> bool {
        other
            .atoms
            .iter()
            .all(|(atom, &count)| self.atoms.get(atom).copied().unwrap_or(0) >= count)
    }

    /// Checked algebra retains negative counts and removes zero entries.
    pub fn checked_add(&self, other: &Self) -> Result<Self> {
        self.combine(other, 1)
    }

    pub fn checked_sub(&self, other: &Self) -> Result<Self> {
        self.combine(other, -1)
    }

    fn combine(&self, other: &Self, direction: i64) -> Result<Self> {
        let mut result = self.clone();
        result.charge = checked_i32(i64::from(self.charge) + direction * i64::from(other.charge))?;
        for (&atom, &count) in &other.atoms {
            let old = result.atoms.get(&atom).copied().unwrap_or(0);
            let next = checked_i32(i64::from(old) + direction * i64::from(count))?;
            if next == 0 {
                result.atoms.remove(&atom);
            } else {
                result.atoms.insert(atom, next);
            }
        }
        Ok(result)
    }

    pub fn checked_scale(&self, factor: i32) -> Result<Self> {
        if factor == 0 {
            return Ok(Self::default());
        }
        let mut result = self.clone();
        result.charge = checked_i32(i64::from(self.charge) * i64::from(factor))?;
        for count in result.atoms.values_mut() {
            *count = checked_i32(i64::from(*count) * i64::from(factor))?;
        }
        Ok(result)
    }
}

fn checked_i32(value: i64) -> Result<i32> {
    i32::try_from(value)
        .map_err(|_| Error::InvalidValue("formula count or charge overflows i32".into()))
}

fn parse_atom(bytes: &[u8], position: &mut usize) -> Result<Atom> {
    let isotope = if bytes.get(*position) == Some(&b'(') {
        *position += 1;
        let start = *position;
        while bytes.get(*position).is_some_and(u8::is_ascii_digit) {
            *position += 1;
        }
        if *position == start || bytes.get(*position) != Some(&b')') {
            return Err(Error::InvalidValue(
                "expected isotope label such as (13)C".into(),
            ));
        }
        let number = std::str::from_utf8(&bytes[start..*position])
            .expect("ASCII digits")
            .parse::<u16>()
            .map_err(|_| Error::InvalidValue("isotope number is too large".into()))?;
        *position += 1;
        Some(number)
    } else {
        None
    };
    let start = *position;
    if !bytes.get(*position).is_some_and(u8::is_ascii_uppercase) {
        return Err(Error::InvalidValue(format!(
            "expected element at byte {start}"
        )));
    }
    *position += 1;
    while bytes.get(*position).is_some_and(u8::is_ascii_lowercase) {
        *position += 1;
    }
    let symbol = std::str::from_utf8(&bytes[start..*position]).expect("ASCII element symbol");
    Atom::resolve(symbol, isotope)
}

impl FromStr for EmpiricalFormula {
    type Err = Error;

    fn from_str(input: &str) -> Result<Self> {
        let input = input.trim();
        let bytes = input.as_bytes();
        let mut position = 0;
        let mut result = Self::default();
        while position < bytes.len() {
            // A sign encountered between complete atom tokens is a terminal charge.
            if matches!(bytes[position], b'+' | b'-') {
                let suffix = &input[position..];
                result.charge = if suffix == "+" {
                    1
                } else if suffix == "-" {
                    -1
                } else {
                    if !suffix.as_bytes()[1..].iter().all(u8::is_ascii_digit) {
                        return Err(Error::InvalidValue("invalid terminal charge".into()));
                    }
                    suffix
                        .parse()
                        .map_err(|_| Error::InvalidValue("invalid terminal charge".into()))?
                };
                return Ok(result);
            }
            let atom = parse_atom(bytes, &mut position)?;
            let count_start = position;
            // A minus before digits directly after a symbol means negative count.
            if bytes.get(position) == Some(&b'-')
                && bytes.get(position + 1).is_some_and(u8::is_ascii_digit)
            {
                position += 1;
            }
            while bytes.get(position).is_some_and(u8::is_ascii_digit) {
                position += 1;
            }
            let count = if position == count_start {
                1
            } else {
                input[count_start..position]
                    .parse::<i32>()
                    .map_err(|_| Error::InvalidValue("atom count overflows i32".into()))?
            };
            let old = result.atoms.get(&atom).copied().unwrap_or(0);
            let next = checked_i32(i64::from(old) + i64::from(count))?;
            if next == 0 {
                result.atoms.remove(&atom);
            } else {
                result.atoms.insert(atom, next);
            }
        }
        Ok(result)
    }
}

impl fmt::Display for EmpiricalFormula {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (atom, count) in &self.atoms {
            if let Some(isotope) = atom.isotope {
                write!(f, "({isotope})")?;
            }
            write!(f, "{}{count}", atom.symbol)?;
        }
        if self.charge != 0 {
            write!(f, "{:+}", self.charge)?;
        }
        Ok(())
    }
}

// Internal residue formulas from ResidueDB.cpp, i.e. full amino acid minus H2O.
// Columns: C, H, N, O, S, Se. J shares the composition of I and L.
fn residue_composition(residue: u8) -> Option<[i32; 6]> {
    Some(match residue {
        b'A' => [3, 5, 1, 1, 0, 0],
        b'C' => [3, 5, 1, 1, 1, 0],
        b'D' => [4, 5, 1, 3, 0, 0],
        b'E' => [5, 7, 1, 3, 0, 0],
        b'F' => [9, 9, 1, 1, 0, 0],
        b'G' => [2, 3, 1, 1, 0, 0],
        b'H' => [6, 7, 3, 1, 0, 0],
        b'I' | b'L' | b'J' => [6, 11, 1, 1, 0, 0],
        b'K' => [6, 12, 2, 1, 0, 0],
        b'M' => [5, 9, 1, 1, 1, 0],
        b'N' => [4, 6, 2, 2, 0, 0],
        b'P' => [5, 7, 1, 1, 0, 0],
        b'Q' => [5, 8, 2, 2, 0, 0],
        b'R' => [6, 12, 4, 1, 0, 0],
        b'S' => [3, 5, 1, 2, 0, 0],
        b'T' => [4, 7, 1, 2, 0, 0],
        b'V' => [5, 9, 1, 1, 0, 0],
        b'W' => [11, 10, 2, 1, 0, 0],
        b'Y' => [9, 9, 1, 2, 0, 0],
        b'U' => [3, 5, 1, 1, 0, 1],
        b'O' => [12, 19, 3, 2, 0, 0],
        _ => return None,
    })
}

fn composition_formula(counts: [i32; 6]) -> EmpiricalFormula {
    let atoms = ["C", "H", "N", "O", "S", "Se"]
        .into_iter()
        .zip(counts)
        .filter(|(_, count)| *count != 0)
        .map(|(symbol, count)| {
            (
                Atom {
                    symbol,
                    isotope: None,
                },
                count,
            )
        })
        .collect();
    EmpiricalFormula { atoms, charge: 0 }
}

fn ion_mz(ion_mass: f64, charge: i32) -> Result<f64> {
    if charge == 0 {
        return Err(Error::InvalidValue(
            "m/z is undefined for zero charge".into(),
        ));
    }
    if !ion_mass.is_finite() || ion_mass < 0.0 {
        return Err(Error::InvalidValue(
            "calculated ion mass must be finite and nonnegative".into(),
        ));
    }
    Ok(ion_mass / f64::from(charge.unsigned_abs()))
}

/// Supported fragment series (neutral b = internal residues, neutral y adds H2O).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IonSeries {
    B,
    Y,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FragmentIon {
    pub series: IonSeries,
    /// Number of residues retained in the fragment.
    pub ordinal: usize,
    pub charge: u8,
    pub mz: f64,
}

pub mod proforma;

pub(crate) mod metabo_elements;
