// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The remainder of `CHEMISTRY/EmpiricalFormula.h` on [`EmpiricalFormula`].
//!
//! Parsing, masses, charge, atom counts, containment and the checked algebra
//! live with the type in [`crate::chemistry`]; this module adds the members
//! that were still missing: construction from one element, the
//! [`EmpiricalFormula::hydrogen`] and [`EmpiricalFormula::water`] factories,
//! element queries, the symbol map, iteration, the lightest-isotope weight, the
//! theoretical isotope count, charge adducts, the two source text forms,
//! ordering and the source hash.
//!
//! # API mapping
//!
//! | Source (`EmpiricalFormula.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `EmpiricalFormula()`, copy and move constructors, `operator=` | `Default`, `Clone` |
//! | `explicit EmpiricalFormula(const std::string&)`, `static fromString` | [`EmpiricalFormula::parse`], `FromStr` |
//! | `EmpiricalFormula(SignedSize number, const Element*, SignedSize charge = 0)` | [`EmpiricalFormula::from_element`] |
//! | `getMonoWeight`, `getAverageWeight` | [`EmpiricalFormula::mono_mass`], [`EmpiricalFormula::average_mass`] |
//! | `getLightestIsotopeWeight` | [`EmpiricalFormula::lightest_isotope_mass`] |
//! | `calculateTheoreticalIsotopesNumber` | [`EmpiricalFormula::theoretical_isotopes_number`] |
//! | `estimateFromWeightAndComp`, `estimateFromMonoWeightAndComp`, `estimateFromWeightAndCompAndS` | [`crate::chemistry::isotopes::AveragineComposition`] |
//! | `getIsotopeDistribution(const IsotopePatternGenerator&)` | the generators' `run` ([`crate::chemistry::isotopes::CoarseIsotopePatternGenerator::run`], fine isotope generator) |
//! | `getConditionalFragmentIsotopeDist` | [`crate::chemistry::isotopes::CoarseIsotopePatternGenerator::calc_fragment_isotope_dist`] |
//! | `getNumberOf(const Element*)` | [`EmpiricalFormula::number_of`]; by symbol [`EmpiricalFormula::count`] |
//! | `getNumberOfAtoms` | [`EmpiricalFormula::atom_count`] |
//! | `getCharge`, `setCharge` | [`EmpiricalFormula::charge`], [`EmpiricalFormula::set_charge`], [`EmpiricalFormula::with_charge`] |
//! | `addChargeAdduct(Int count = 1, const std::string& adduct = "H")` | [`EmpiricalFormula::add_charge_adduct`] |
//! | `toString` | [`EmpiricalFormula::to_source_string`] |
//! | `toMap` | [`EmpiricalFormula::to_map`] |
//! | `operator+=`, `operator+`, `operator-=`, `operator-`, `operator*` | [`EmpiricalFormula::checked_add`], [`EmpiricalFormula::checked_sub`], [`EmpiricalFormula::checked_scale`] |
//! | `isEmpty`, `isCharged`, `hasElement(const Element*)`, `contains` | [`EmpiricalFormula::is_empty`], [`EmpiricalFormula::is_charged`], [`EmpiricalFormula::has_element`], [`EmpiricalFormula::contains`] |
//! | `operator==`, `operator!=` | `PartialEq`, `Eq` |
//! | `operator<` | `Ord` |
//! | `friend operator<<(std::ostream&, const EmpiricalFormula&)` | [`EmpiricalFormula::source_stream_text`]; `Display` is the native round-trip form |
//! | `begin()`, `end()` (const and mutable), the four iterator typedefs | [`EmpiricalFormula::entries`] |
//! | `static hydrogen(int n_atoms = 1)`, `static water(int n_molecules = 1)` | [`EmpiricalFormula::hydrogen`], [`EmpiricalFormula::water`] |
//! | `std::hash<EmpiricalFormula>` | `Hash`; the source's value is [`EmpiricalFormula::source_hash64`] |
//! | protected `removeZeroedElements_`, `parseFormula_`, `formula_`, `charge_` | every operation keeps zero counts out; the parser; private fields |
//!
//! # Preserved source conventions
//!
//! - An isotope-labelled atom is its own entry, as the source's `(13)C`
//!   `ElementDB` element is its own map key: natural carbon and `(13)C` are
//!   counted, queried and hashed separately. `D` and `T` are `(2)H` and `(3)H`.
//! - [`EmpiricalFormula::to_map`], [`EmpiricalFormula::to_source_string`] and
//!   [`EmpiricalFormula::source_stream_text`] key by the source symbol text
//!   (`(13)C`, `C`, `H`) in byte order, as the source's `std::map<std::string, …>`
//!   does, so `(13)C` sorts before `C`.
//! - [`EmpiricalFormula::add_charge_adduct`] ignores any charge written in the
//!   adduct, scales it by a signed count, adds it and sets the charge to zero.
//! - [`EmpiricalFormula::source_hash64`] is the source's hash value: entries
//!   sorted by symbol, each symbol's FNV-1a hash and its count's `hash_int` over
//!   the eight bytes of `SignedSize`, then the charge's `hash_int` over the four
//!   bytes of `Int`.
//!
//! # Native differences
//!
//! - The source orders with `operator<` on its `std::map<const Element*, …>`,
//!   whose third criterion compares entries in element *address* order, which
//!   varies between runs. `Ord` keeps the source's first two criteria (number
//!   of entries, then charge) and compares entries in the fixed order of this
//!   type (symbol, then natural before labelled, then mass number).
//! - [`EmpiricalFormula::from_element`] with a zero count stores nothing; the
//!   source stores a zero entry, so its `isEmpty` is false for that formula.
//! - [`EmpiricalFormula::theoretical_isotopes_number`] refuses a negative atom
//!   count on the binomial branch, which the source converts to an unsigned
//!   32-bit value, and returns an exact binomial coefficient correctly rounded
//!   to `f64` where Boost evaluates a beta function above 170 atoms (identical
//!   below `2^53`).
//! - Counts are `i32` and arithmetic is checked: [`EmpiricalFormula::water`]
//!   and [`EmpiricalFormula::add_charge_adduct`] return an error instead of
//!   overflowing.

use super::{Atom, Element, EmpiricalFormula, element, element_table, parse_atom};
use crate::concept::hash_utils::{fnv1a_hash_string, hash_combine, hash_int};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;

/// One entry of an [`EmpiricalFormula`]: an element, or one of its labelled
/// isotopes, with its signed count.
///
/// Source: one `(const Element*, SignedSize)` pair of the formula's map, as its
/// `begin()`/`end()` iterators yield them. A labelled isotope is a separate
/// `ElementDB` element in the source; here it is the natural element with the
/// isotope's mass number.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FormulaEntry {
    /// The natural element.
    pub element: &'static Element,
    /// The mass number of a labelled isotope (`13` for `(13)C`), or `None` for
    /// the natural element.
    pub isotope: Option<u16>,
    /// The signed number of atoms; never zero.
    pub count: i32,
}

impl FormulaEntry {
    /// The source `ElementDB` symbol of this entry: `C` for natural carbon,
    /// `(13)C` for carbon-13 (also for a formula written with `D` or `T`, which
    /// are `(2)H` and `(3)H`).
    pub fn symbol(&self) -> String {
        source_symbol(self.element.symbol(), self.isotope)
    }
}

impl EmpiricalFormula {
    /// `count` atoms of one element with the given charge: source
    /// `EmpiricalFormula(SignedSize number, const Element* element, SignedSize charge)`.
    ///
    /// An isotope element of [`super::element_db`] (`(13)C`) gives the labelled
    /// isotope, as the source's isotope element is its own map key. A zero count leaves the formula without atoms (see the module notes).
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::{element, EmpiricalFormula};
    ///
    /// let carbon = element("C").unwrap();
    /// let c6 = EmpiricalFormula::from_element(6, carbon, 1);
    /// assert_eq!(c6, EmpiricalFormula::parse("C6+").unwrap());
    /// ```
    pub fn from_element(count: i32, element: &Element, charge: i32) -> Self {
        let mut atoms = BTreeMap::new();
        if count != 0 {
            atoms.insert(atom_of(element), count);
        }
        Self { atoms, charge }
    }

    /// `n_atoms` hydrogen atoms, uncharged: source `static hydrogen(int n_atoms = 1)`.
    pub fn hydrogen(n_atoms: i32) -> Self {
        Self::from_element(n_atoms, hydrogen(), 0)
    }

    /// `n_molecules` water molecules (`H2O` each), uncharged: source
    /// `static water(int n_molecules = 1)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `2 * n_molecules` overflows `i32`;
    /// the source overflows a signed `int`.
    pub fn water(n_molecules: i32) -> Result<Self> {
        let hydrogens = n_molecules
            .checked_mul(2)
            .ok_or_else(|| Error::InvalidValue("water: hydrogen count overflows i32".into()))?;
        let mut formula = Self::from_element(hydrogens, hydrogen(), 0);
        if n_molecules != 0 {
            formula.atoms.insert(
                Atom {
                    symbol: "O",
                    isotope: None,
                },
                n_molecules,
            );
        }
        Ok(formula)
    }

    /// Replaces the charge in place: source `setCharge(Int)`.
    pub fn set_charge(&mut self, charge: i32) {
        self.charge = charge;
    }

    /// Whether the charge is nonzero: source `isCharged`.
    pub fn is_charged(&self) -> bool {
        self.charge != 0
    }

    /// Whether the formula has an entry for this element: source
    /// `hasElement(const Element*)`. A natural element does not match its
    /// labelled isotopes, as the source's `(13)C` element is a different
    /// pointer; pass the isotope element of [`super::element_db`] for those.
    pub fn has_element(&self, element: &Element) -> bool {
        self.number_of(element) != 0
    }

    /// The count of this element, zero when absent: source
    /// `getNumberOf(const Element*)`. A natural element does not count its
    /// labelled isotopes; pass the isotope element of [`super::element_db`]
    /// (`(13)C`), or use [`Self::count`] with the symbol.
    pub fn number_of(&self, element: &Element) -> i32 {
        self.atoms.get(&atom_of(element)).copied().unwrap_or(0)
    }

    /// Every entry in the order of this type (symbol, natural before labelled,
    /// ascending mass number): source `begin()`/`end()`.
    ///
    /// The source iterates its map in element address order, which is not
    /// fixed. Entries are read-only; build a changed formula with the checked
    /// algebra.
    pub fn entries(&self) -> impl ExactSizeIterator<Item = FormulaEntry> + '_ {
        self.atoms.iter().map(|(atom, &count)| FormulaEntry {
            element: atom_element(*atom),
            isotope: atom.isotope,
            count,
        })
    }

    /// Source symbol to count, in symbol byte order: source `toMap`.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::EmpiricalFormula;
    ///
    /// let map = EmpiricalFormula::parse("C2H5(13)C").unwrap().to_map();
    /// let pairs: Vec<_> = map.iter().map(|(s, n)| (s.as_str(), *n)).collect();
    /// assert_eq!(pairs, [("(13)C", 1), ("C", 2), ("H", 5)]);
    /// ```
    pub fn to_map(&self) -> BTreeMap<String, i32> {
        self.atoms
            .iter()
            .map(|(atom, &count)| (source_symbol(atom.symbol, atom.isotope), count))
            .collect()
    }

    /// Source `toString`: each symbol of [`Self::to_map`] followed by its count,
    /// written even when it is `1` or negative; the charge is not written.
    ///
    /// `C2H5+` gives `C2H5`; `H-2O` gives `H-2O1`.
    pub fn to_source_string(&self) -> String {
        let mut text = String::new();
        for (symbol, count) in self.to_map() {
            text.push_str(&symbol);
            text.push_str(&count.to_string());
        }
        text
    }

    /// The text source `operator<<(std::ostream&, const EmpiricalFormula&)`
    /// writes.
    ///
    /// Each symbol of [`Self::to_map`] is followed by its count only when the
    /// count exceeds one, so a count of one and every negative count print the
    /// bare symbol. A charge of `+1` or `-1` prints `+` or `-`; any other
    /// positive charge prints `+` and the number; any other negative charge
    /// prints `-` followed by the negative number, so `-2` prints `--2`. That is
    /// the source's behavior, kept for fidelity; [`std::fmt::Display`] writes
    /// the round-trip form instead.
    pub fn source_stream_text(&self) -> String {
        let mut text = String::new();
        for (symbol, count) in self.to_map() {
            text.push_str(&symbol);
            if count > 1 {
                text.push_str(&count.to_string());
            }
        }
        match self.charge {
            0 => {}
            1 => text.push('+'),
            -1 => text.push('-'),
            charge if charge > 0 => {
                text.push('+');
                text.push_str(&charge.to_string());
            }
            charge => {
                text.push('-');
                text.push_str(&charge.to_string());
            }
        }
        text
    }

    /// Mass of the formula built from each element's lightest tabulated isotope
    /// (a labelled isotope's own mass), plus `charge * PROTON_MASS_U`: source
    /// `getLightestIsotopeWeight`.
    pub fn lightest_isotope_mass(&self) -> f64 {
        super::isotopes::lightest_mass(self)
    }

    /// Source `calculateTheoreticalIsotopesNumber`: the product over entries of
    /// the number of isotope combinations the source counts.
    ///
    /// For each entry, `k` is the number of isotopes with a nonzero natural
    /// abundance (one for a labelled isotope). When `k > 1` and the count is not
    /// one, the factor is the binomial coefficient `C(count, k)`; otherwise it
    /// is `count * k`. An empty formula gives `1`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when the binomial branch meets a negative
    /// count (the source wraps it to an unsigned value), a count smaller than
    /// `k` (Boost raises a domain error), or a coefficient beyond the `f64`
    /// range (Boost raises an overflow error).
    pub fn theoretical_isotopes_number(&self) -> Result<f64> {
        let mut total = 1.0f64;
        for (atom, &count) in &self.atoms {
            let non_trace = match atom.isotope {
                Some(_) => 1u32,
                None => {
                    let isotopes = atom_element(*atom).isotopes();
                    // At most a few dozen isotopes per element.
                    isotopes.iter().filter(|i| i.abundance != 0.0).count() as u32
                }
            };
            if non_trace > 1 && count != 1 {
                let n = u32::try_from(count).map_err(|_| {
                    Error::InvalidValue(format!(
                        "theoretical isotope number: negative count {count} for {}",
                        source_symbol(atom.symbol, atom.isotope)
                    ))
                })?;
                total *= binomial(n, non_trace)?;
            } else {
                total *= (i64::from(count) * i64::from(non_trace)) as f64;
            }
        }
        Ok(total)
    }

    /// Adds `count` copies of `adduct` and makes the formula neutral: source
    /// `addChargeAdduct(Int count = 1, const std::string& adduct = "H")`.
    ///
    /// A charge written in `adduct` is ignored; a negative `count` removes
    /// atoms. The charge is set to zero unconditionally, so to neutralise an
    /// already-charged formula consistently pass its charge as `count`. The
    /// source's defaults are `add_charge_adduct(1, "H")`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidValue`] when `adduct` is not a formula or a count
    /// overflows `i32`; the formula is then unchanged.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::EmpiricalFormula;
    ///
    /// let mut glucose = EmpiricalFormula::parse("C6H12O6").unwrap();
    /// glucose.add_charge_adduct(2, "H").unwrap();
    /// assert_eq!(glucose, EmpiricalFormula::parse("C6H14O6").unwrap());
    /// assert_eq!(glucose.charge(), 0);
    /// ```
    pub fn add_charge_adduct(&mut self, count: i32, adduct: &str) -> Result<&mut Self> {
        let adduct = EmpiricalFormula::parse(adduct)?.with_charge(0);
        let mut result = self.checked_add(&adduct.checked_scale(count)?)?;
        result.charge = 0;
        *self = result;
        Ok(self)
    }

    /// The value source `std::hash<EmpiricalFormula>` computes.
    ///
    /// Entries are sorted by source symbol; for each, the symbol's
    /// [`fnv1a_hash_string`] and the count's [`hash_int`] over the eight bytes
    /// of `SignedSize` are combined into a seed starting at zero, then the
    /// charge's `hash_int` over the four bytes of `Int`.
    pub fn source_hash64(&self) -> u64 {
        let mut seed = 0u64;
        for (symbol, count) in self.to_map() {
            hash_combine(&mut seed, fnv1a_hash_string(&symbol));
            hash_combine(&mut seed, hash_int(i64::from(count)));
        }
        hash_combine(&mut seed, hash_int(self.charge));
        seed
    }
}

impl PartialOrd for EmpiricalFormula {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for EmpiricalFormula {
    /// Source `operator<`: fewer entries first, then the smaller charge, then
    /// the entries compared in order (see the module notes for the order).
    fn cmp(&self, other: &Self) -> Ordering {
        self.graph_cmp(other)
    }
}

/// The formula key of a database element: the natural element, or for an
/// isotope element of [`super::element_db`] (`(13)C`) the labelled isotope.
fn atom_of(element: &Element) -> Atom {
    let mut position = 0;
    // Elements exist only as entries of the compiled table and the isotope
    // elements derived from it, whose symbols always parse.
    parse_atom(element.symbol().as_bytes(), &mut position).unwrap_or(Atom {
        symbol: element.symbol(),
        isotope: None,
    })
}

fn hydrogen() -> &'static Element {
    &element_table()[0]
}

fn atom_element(atom: Atom) -> &'static Element {
    // Atoms are only constructed from entries of the static element table.
    element(atom.symbol).expect("validated formula element")
}

fn source_symbol(symbol: &str, isotope: Option<u16>) -> String {
    match isotope {
        Some(number) => format!("({number}){symbol}"),
        None => symbol.to_owned(),
    }
}

/// Exact `C(n, k)` rounded to `f64`, falling back to a floating product once
/// the exact value leaves `u128`.
fn binomial(n: u32, k: u32) -> Result<f64> {
    if k > n {
        return Err(Error::InvalidValue(format!(
            "theoretical isotope number: C({n}, {k}) is undefined for k > n"
        )));
    }
    let k = k.min(n - k);
    let mut exact: Option<u128> = Some(1);
    let mut approximate = 1.0f64;
    for i in 0..k {
        let numerator = u128::from(n - i);
        let denominator = u128::from(i + 1);
        // The running product C(n, i + 1) is exactly divisible at every step.
        exact = exact
            .and_then(|value| value.checked_mul(numerator))
            .map(|v| v / denominator);
        approximate = approximate * f64::from(n - i) / f64::from(i + 1);
    }
    let value = match exact {
        Some(value) => value as f64,
        None => approximate.round(),
    };
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidValue(format!(
            "theoretical isotope number: C({n}, {k}) overflows f64"
        )))
    }
}
