// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The element database: `CHEMISTRY/ElementDB.h`, and the remaining
//! value-type members of `CHEMISTRY/Element.h` on [`Element`].
//!
//! The source `ElementDB` is a singleton holding every element the source
//! declares, reachable by name, by symbol and by atomic number. Besides the 84
//! natural elements of [`element_table`], it holds one element per tabulated
//! isotope, named and symbolised with the mass number in brackets
//! (`(13)Carbon`, `(13)C`), whose only isotope is that one at abundance one and
//! whose average and monoisotopic weights are both its mass. `D` and `T` are
//! additional symbols for `(2)H` and `(3)H`. Here the database is built once,
//! on first use, and never changes; every element it returns lives for the
//! whole program, so two lookups of one element return the same reference, as
//! the source returns the same pointer.
//!
//! # API mapping
//!
//! | Source (`ElementDB.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `static getInstance()` | the free functions of this module |
//! | `getNames()` | [`names`](crate::chemistry::element_db::names) |
//! | `getSymbols()` | [`symbols`](crate::chemistry::element_db::symbols) |
//! | `getAtomicNumbers()` | [`atomic_numbers`](crate::chemistry::element_db::atomic_numbers) |
//! | `getElement(const std::string&)` | [`get_element`](crate::chemistry::element_db::get_element) |
//! | `getElement(unsigned int)` | [`get_element_by_atomic_number`](crate::chemistry::element_db::get_element_by_atomic_number) |
//! | `hasElement(const std::string&)` | [`has_element`](crate::chemistry::element_db::has_element) |
//! | `hasElement(unsigned int)` | [`has_atomic_number`](crate::chemistry::element_db::has_atomic_number) |
//! | protected `parseIsotopeDistribution_`, `calculateAvgWeight_`, `calculateMonoWeight_`, `storeElements_`, `buildElement_`, `addElementToMaps_`, `storeIsotopes_`, `clear_`; the private constructors | the compiled table and the one-time build of this module |
//!
//! | Source (`Element.h` at core `bc9cc12`) | Rust |
//! |---|---|
//! | `getIsotopeDistribution()` | [`Element::isotope_distribution`]; the table as [`Element::isotopes`] |
//! | `operator<` | `Ord` |
//! | `friend operator<<(std::ostream&, const Element&)` | [`Element::source_stream_text`] |
//! | `std::hash<Element>` | `Hash`; the source's value is [`Element::source_hash64`] |
//!
//! # Differences from the source
//!
//! - Iridium carries its declared iridium isotopes, where the pinned source
//!   builds it from rhenium's tables (CPP-249); so its isotope elements are
//!   `(191)Ir` and `(193)Ir`, where the source has `(185)Ir` and `(187)Ir`.
//! - [`super::element`](crate::chemistry::element) remains the lookup of natural elements only; use
//!   [`get_element`](crate::chemistry::element_db::get_element) for the source's lookup, which also finds the isotope
//!   elements and the `D`/`T` aliases.

use super::{Element, Isotope, element_table};
use crate::chemistry::isotopes::{IsotopeDistribution, IsotopePeak};
use crate::concept::fuzzy_string_comparator::format_g;
use crate::concept::hash_utils::{fnv1a_hash_string, hash_combine, hash_float, hash_int};
use crate::concept::numeric_formatting::append_numeric_f64;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::sync::OnceLock;

struct Maps {
    names: HashMap<&'static str, &'static Element>,
    symbols: HashMap<&'static str, &'static Element>,
    atomic_numbers: HashMap<u32, &'static Element>,
}

/// The isotope elements, built once. Their text and one-entry isotope tables
/// are leaked on that one build so that they are `'static` like the compiled
/// natural table; the total is a few hundred short strings.
fn isotope_elements() -> &'static [Element] {
    static ISOTOPE_ELEMENTS: OnceLock<Vec<Element>> = OnceLock::new();
    ISOTOPE_ELEMENTS.get_or_init(|| {
        let mut elements = Vec::new();
        for natural in element_table() {
            for isotope in natural.isotopes() {
                let prefix = format!("({})", isotope.mass_number);
                let isotopes: &'static [Isotope] = Box::leak(Box::new([Isotope {
                    mass_number: isotope.mass_number,
                    mass: isotope.mass,
                    abundance: 1.0,
                }]));
                elements.push(Element {
                    name: Box::leak(format!("{prefix}{}", natural.name()).into_boxed_str()),
                    symbol: Box::leak(format!("{prefix}{}", natural.symbol()).into_boxed_str()),
                    atomic_number: natural.atomic_number(),
                    isotopes,
                });
            }
        }
        elements
    })
}

fn maps() -> &'static Maps {
    static MAPS: OnceLock<Maps> = OnceLock::new();
    MAPS.get_or_init(|| {
        let mut maps = Maps {
            names: HashMap::new(),
            symbols: HashMap::new(),
            atomic_numbers: HashMap::new(),
        };
        for natural in element_table() {
            maps.names.insert(natural.name(), natural);
            maps.symbols.insert(natural.symbol(), natural);
            maps.atomic_numbers
                .insert(u32::from(natural.atomic_number()), natural);
        }
        for isotope in isotope_elements() {
            maps.names.insert(isotope.name(), isotope);
            maps.symbols.insert(isotope.symbol(), isotope);
        }
        for (alias, target) in [("D", "(2)H"), ("T", "(3)H")] {
            if let Some(&found) = maps.symbols.get(target) {
                maps.symbols.insert(alias, found);
            }
        }
        maps
    })
}

/// Every element by name, natural and isotope (`Carbon`, `(13)Carbon`):
/// source `getNames()`.
pub fn names() -> &'static HashMap<&'static str, &'static Element> {
    &maps().names
}

/// Every element by symbol, natural and isotope (`C`, `(13)C`), plus `D` and
/// `T` for `(2)H` and `(3)H`: source `getSymbols()`.
pub fn symbols() -> &'static HashMap<&'static str, &'static Element> {
    &maps().symbols
}

/// The natural elements by atomic number: source `getAtomicNumbers()`. Isotope
/// elements are not in this map, as in the source.
pub fn atomic_numbers() -> &'static HashMap<u32, &'static Element> {
    &maps().atomic_numbers
}

/// The element with this symbol, or else this name: source
/// `getElement(const std::string&)`, which returns `nullptr` where this returns
/// `None`.
///
/// # Examples
///
/// ```
/// use openms::chemistry::element_db::get_element;
///
/// let hydrogen = get_element("Hydrogen").unwrap();
/// assert!(std::ptr::eq(hydrogen, get_element("H").unwrap()));
/// let deuterium = get_element("D").unwrap();
/// assert_eq!(deuterium.name(), "(2)Hydrogen");
/// assert_eq!(deuterium.symbol(), "(2)H");
/// assert!(get_element("Xx").is_none());
/// ```
pub fn get_element(name: &str) -> Option<&'static Element> {
    let maps = maps();
    maps.symbols
        .get(name)
        .or_else(|| maps.names.get(name))
        .copied()
}

/// The natural element with this atomic number: source
/// `getElement(unsigned int)`, which returns `nullptr` where this returns
/// `None`.
pub fn get_element_by_atomic_number(atomic_number: u32) -> Option<&'static Element> {
    maps().atomic_numbers.get(&atomic_number).copied()
}

/// Whether an element has this name or symbol: source
/// `hasElement(const std::string&)`.
pub fn has_element(name: &str) -> bool {
    get_element(name).is_some()
}

/// Whether a natural element has this atomic number: source
/// `hasElement(unsigned int)`.
pub fn has_atomic_number(atomic_number: u32) -> bool {
    maps().atomic_numbers.contains_key(&atomic_number)
}

impl Element {
    /// The isotope table as an [`IsotopeDistribution`] of (mass, abundance)
    /// peaks in mass-number order: source `getIsotopeDistribution()`.
    ///
    /// The source stores each abundance in `Peak1D`'s `float` intensity; this
    /// keeps the tabulated `f64`.
    pub fn isotope_distribution(&self) -> IsotopeDistribution {
        IsotopeDistribution::from_peaks(
            self.isotopes()
                .iter()
                .map(|isotope| IsotopePeak {
                    mass: isotope.mass,
                    probability: isotope.abundance,
                })
                .collect(),
        )
        // Tabulated masses and abundances are finite and nonnegative, and no
        // element has more isotopes than the distribution limit.
        .unwrap_or_else(|_| IsotopeDistribution::empty())
    }

    /// Whether this is one of the database's isotope elements, such as
    /// `(13)C`, rather than a natural element.
    pub fn is_isotope(&self) -> bool {
        self.symbol().starts_with('(')
    }

    /// The text source `operator<<(std::ostream&, const Element&)` writes.
    ///
    /// Name, symbol, atomic number, average weight and monoisotopic weight,
    /// separated by single spaces, the two weights in the stream's default
    /// six significant digits (`%g`); then, for every isotope whose `float`
    /// abundance is positive, a space, its mass in full precision (the
    /// `DPosition` output, [`crate::concept::precision_wrapper`]), `=`, the
    /// `float` abundance times 100 in six significant digits, and `%`.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::element;
    ///
    /// assert_eq!(
    ///     element("C").unwrap().source_stream_text(),
    ///     "Carbon C 6 12.0107 12 12.0=98.93% 13.003355000000001=1.07%"
    /// );
    /// ```
    pub fn source_stream_text(&self) -> String {
        let mut text = format!(
            "{} {} {} {} {}",
            self.name(),
            self.symbol(),
            self.atomic_number(),
            format_g(self.average_mass()),
            format_g(self.mono_mass())
        );
        for isotope in self.isotopes() {
            let abundance = isotope.abundance as f32;
            if abundance > 0.0 {
                text.push(' ');
                append_numeric_f64(isotope.mass, &mut text, 15, false);
                text.push('=');
                text.push_str(&format_g(f64::from(abundance * 100.0)));
                text.push('%');
            }
        }
        text
    }

    /// The value source `std::hash<Element>` computes.
    ///
    /// The seed is the [`fnv1a_hash_string`] of the name; combined into it, in
    /// order, are the symbol's, `hash_int` of the atomic number as a 4-byte
    /// `unsigned int`, `hash_float` of the average and of the monoisotopic
    /// weight, and the [`IsotopeDistribution::source_hash64`] of the isotope
    /// distribution (abundances narrowed to `float` as the source stores them).
    pub fn source_hash64(&self) -> u64 {
        let mut seed = fnv1a_hash_string(self.name());
        hash_combine(&mut seed, fnv1a_hash_string(self.symbol()));
        hash_combine(&mut seed, hash_int(u32::from(self.atomic_number())));
        hash_combine(&mut seed, hash_float(self.average_mass()));
        hash_combine(&mut seed, hash_float(self.mono_mass()));
        hash_combine(&mut seed, self.isotope_distribution().source_hash64());
        seed
    }
}

impl Eq for Element {}

impl PartialOrd for Element {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Element {
    /// Source `operator<`: atomic number, monoisotopic weight, symbol, name,
    /// average weight, then the isotope distributions.
    fn cmp(&self, other: &Self) -> Ordering {
        self.atomic_number()
            .cmp(&other.atomic_number())
            .then_with(|| compare(self.mono_mass(), other.mono_mass()))
            .then_with(|| self.symbol().cmp(other.symbol()))
            .then_with(|| self.name().cmp(other.name()))
            .then_with(|| compare(self.average_mass(), other.average_mass()))
            .then_with(|| {
                self.isotope_distribution()
                    .cmp(&other.isotope_distribution())
            })
    }
}

impl Hash for Element {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name().hash(state);
        self.symbol().hash(state);
        self.atomic_number().hash(state);
        self.isotope_distribution().hash(state);
    }
}

fn compare(left: f64, right: f64) -> Ordering {
    left.partial_cmp(&right).unwrap_or(Ordering::Equal)
}
