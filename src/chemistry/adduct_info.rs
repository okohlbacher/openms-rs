// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Ionizing adduct formulas and electron-aware mass arithmetic, following
//! pinned OpenMS `7c029e8` AdductInfo. The adduct formula stays neutral; the
//! signed ion charge and molecular multiplier are stored separately.

use super::{ELECTRON_MASS_U, EmpiricalFormula, PROTON_MASS_U};
use crate::{Error, Result};
use std::str::FromStr;

pub const MAX_ADDUCT_TEXT_BYTES: usize = 65_536;
pub const MAX_ADDUCT_TERMS: usize = 4_096;
pub const MAX_ADDUCT_WORK: usize = 1_000_000;

/// A named neutral atomic adduct, signed ion charge and positive n-mer factor.
/// Names participate in equality: equivalent chemical expressions need not be
/// equal records. Fields are immutable after checked construction.
#[derive(Clone, Debug)]
pub struct AdductInfo {
    name: String,
    formula: EmpiricalFormula,
    mono_mass: f64,
    charge: i32,
    mol_multiplier: u32,
}

impl PartialEq for AdductInfo {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.formula == other.formula
            && self.charge == other.charge
            && self.mol_multiplier == other.mol_multiplier
    }
}
impl Eq for AdductInfo {}

impl AdductInfo {
    /// Construct from components. The name may be any bounded text, including
    /// empty text; it need not itself be parseable adduct notation. The source
    /// constructor's default molecular multiplier is one.
    pub fn new(
        name: impl AsRef<str>,
        adduct: EmpiricalFormula,
        charge: i32,
        mol_multiplier: u32,
    ) -> Result<Self> {
        let name = name.as_ref();
        check_text(name)?;
        if charge == 0 || charge == i32::MIN {
            return Err(invalid(
                "adduct charge must have a nonzero representable i32 magnitude",
            ));
        }
        if adduct.charge() != 0 {
            return Err(invalid("adduct empirical formula must have zero charge"));
        }
        if mol_multiplier == 0 {
            return Err(invalid("adduct molecular multiplier must be positive"));
        }
        let mono_mass = finite(adduct.mono_mass(), "adduct monoisotopic mass")?;
        Ok(Self {
            name: name.to_owned(),
            formula: adduct,
            mono_mass,
            charge,
            mol_multiplier,
        })
    }

    /// Parse `M+H;1+`, `2M+Na;1+`, etc. Only source SP/TAB/LF/CR whitespace
    /// is removed. The stripped spelling is retained as the record name.
    ///
    /// The final charge sign overrides a signed magnitude. Every formula-side
    /// plus/minus separates terms: `M+H-2;1+` adds H and an empty numeric-only
    /// term, rather than removing two H atoms. Zero coefficients still require
    /// a valid following formula. Bracketed `[M+H]+` notation is not accepted.
    pub fn parse(input: &str) -> Result<Self> {
        check_text(input)?;
        let mut work = MAX_ADDUCT_WORK;
        consume(&mut work, input.len().saturating_mul(4))?;
        let name: String = input
            .chars()
            .filter(|c| !matches!(c, ' ' | '\t' | '\n' | '\r'))
            .collect();
        let (molecular, charge_text) = name
            .split_once(';')
            .filter(|(_, charge)| !charge.contains(';'))
            .ok_or_else(|| invalid("adduct requires exactly one formula/charge semicolon"))?;
        let (polarity, magnitude) = charge_text
            .as_bytes()
            .split_last()
            .filter(|(last, _)| matches!(last, b'+' | b'-'))
            .ok_or_else(|| invalid("adduct charge requires a final plus or minus sign"))?;
        // Charge text is still borrowed from a valid string; +/- is one byte.
        let magnitude = parse_source_i32(&charge_text[..magnitude.len()])?;
        let magnitude = magnitude
            .checked_abs()
            .ok_or_else(|| invalid("adduct charge magnitude overflows i32"))?;
        let charge = if *polarity == b'-' {
            -magnitude
        } else {
            magnitude
        };

        let bytes = molecular.as_bytes();
        let sign = |byte: u8| matches!(byte, b'+' | b'-');
        if bytes.first().is_some_and(|b| sign(*b))
            || bytes.last().is_some_and(|b| sign(*b))
            || bytes.windows(2).any(|pair| sign(pair[0]) && sign(pair[1]))
            || bytes.contains(&b'%')
        {
            return Err(invalid(
                "adduct operators require intervening formulas; percent is forbidden",
            ));
        }
        let first_end = bytes.iter().position(|b| sign(*b)).unwrap_or(bytes.len());
        let molecule = molecular[..first_end]
            .strip_suffix('M')
            .ok_or_else(|| invalid("adduct first term must be M with an optional multiplier"))?;
        let multiplier = if molecule.is_empty() {
            1
        } else {
            parse_source_i32(molecule)?
        };
        let multiplier = u32::try_from(multiplier)
            .map_err(|_| invalid("adduct molecular multiplier must be positive"))?;
        let mut formula = EmpiricalFormula::default();
        let mut position = first_end;
        let mut terms = 0;
        let mut atom_token_bound = 0usize;
        while position < bytes.len() {
            terms += 1;
            if terms > MAX_ADDUCT_TERMS {
                return Err(invalid("adduct term limit exceeded"));
            }
            let plus = bytes[position] == b'+';
            position += 1;
            let start = position;
            while position < bytes.len() && !sign(bytes[position]) {
                position += 1;
            }
            let term = &molecular[start..position];
            let digits = term.bytes().take_while(u8::is_ascii_digit).count();
            let factor = if digits == 0 {
                1
            } else {
                parse_source_i32(&term[..digits])?
            };
            let fragment = &term[digits..];
            // Prevent EmpiricalFormula's broader outer trim from accepting
            // whitespace which the source adduct grammar does not remove.
            if fragment.chars().any(char::is_whitespace) {
                return Err(invalid("unsupported whitespace in adduct formula"));
            }
            let atoms = fragment.bytes().filter(u8::is_ascii_uppercase).count();
            atom_token_bound = atom_token_bound
                .checked_add(atoms)
                .ok_or_else(|| invalid("adduct work accounting overflow"))?;
            // Upper-bound surviving formula keys by parsed atom tokens. This
            // covers repeated checked formula clones/combines without depending
            // on private EmpiricalFormula storage or allocating a token graph.
            let levels = usize::BITS as usize - atom_token_bound.leading_zeros() as usize + 1;
            let cost = fragment
                .len()
                .checked_add(atom_token_bound)
                .and_then(|count| count.checked_mul(8 * levels))
                .and_then(|count| count.checked_add(1))
                .ok_or_else(|| invalid("adduct work accounting overflow"))?;
            consume(&mut work, cost)?;
            let part = EmpiricalFormula::parse(fragment)?.checked_scale(factor)?;
            formula = if plus {
                formula.checked_add(&part)?
            } else {
                formula.checked_sub(&part)?
            };
        }
        Self::new(name, formula, charge, multiplier)
    }

    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn empirical_formula(&self) -> &EmpiricalFormula {
        &self.formula
    }
    pub fn charge(&self) -> i32 {
        self.charge
    }
    pub fn mol_multiplier(&self) -> u32 {
        self.mol_multiplier
    }

    /// Neutral monomer mass from observed m/z, retaining source operation order.
    /// Any finite signed input is accepted; source arithmetic does not enforce
    /// its documentation's positive observed-m/z convention.
    pub fn neutral_mass(&self, observed_mz: f64) -> Result<f64> {
        finite(observed_mz, "observed m/z")?;
        let mass = finite(
            observed_mz * f64::from(self.charge.unsigned_abs()),
            "decharged ion mass",
        )?;
        let mass = finite(mass - self.mono_mass, "adduct-subtracted mass")?;
        let mass = finite(
            mass + f64::from(self.charge) * ELECTRON_MASS_U,
            "electron-corrected mass",
        )?;
        finite(
            mass / f64::from(self.mol_multiplier),
            "neutral monomer mass",
        )
    }

    /// Observed m/z from neutral monomer mass. Charge removes/adds electrons;
    /// no automatic proton or hydrogen atom is added to the neutral formula.
    pub fn mz(&self, neutral_mass: f64) -> Result<f64> {
        finite(neutral_mass, "neutral monomer mass")?;
        let mass = finite(
            neutral_mass * f64::from(self.mol_multiplier),
            "neutral n-mer mass",
        )?;
        let mass = finite(mass + self.mono_mass, "adduct ion mass")?;
        let mass = finite(
            mass - f64::from(self.charge) * ELECTRON_MASS_U,
            "electron-corrected mass",
        )?;
        finite(mass / f64::from(self.charge.unsigned_abs()), "observed m/z")
    }

    /// Source proton-compensated shift: adduct atomic mass minus charge times
    /// `(PROTON_MASS_U + ELECTRON_MASS_U)`. Even average mode retains these
    /// constants; neither mode substitutes the tabulated hydrogen atomic mass.
    pub fn mass_shift(&self, use_average_mass: bool) -> Result<f64> {
        let mass = if use_average_mass {
            self.formula.average_mass()
        } else {
            self.mono_mass
        };
        finite(
            mass - f64::from(self.charge) * (PROTON_MASS_U + ELECTRON_MASS_U),
            "adduct mass shift",
        )
    }

    /// Source signed element-count predicate `candidate.contains(-adduct)`.
    /// Candidate charge and n-mer multiplier are ignored. Plus-only adducts
    /// always pass for nonnegative candidates, but signed candidates can fail.
    pub fn is_compatible(&self, candidate: &EmpiricalFormula) -> bool {
        // The only possible scaling failure is an i32::MIN atom count. Its
        // removal requires 2147483648 atoms, beyond every native i32 candidate,
        // so the source comparison is necessarily false in that case.
        self.formula
            .checked_scale(-1)
            .is_ok_and(|needed| candidate.contains(&needed))
    }
}

impl FromStr for AdductInfo {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}

fn parse_source_i32(text: &str) -> Result<i32> {
    // StringUtils strips one '+' before std::from_chars, which itself accepts
    // a leading '-'. Thus '+-2' is accepted, but '++2' is not.
    let number = text.strip_prefix('+').unwrap_or(text);
    let digits = number.strip_prefix('-').unwrap_or(number);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(invalid("invalid adduct integer"));
    }
    number
        .parse::<i32>()
        .map_err(|_| invalid("adduct integer overflows i32"))
}
fn check_text(text: &str) -> Result<()> {
    if text.len() > MAX_ADDUCT_TEXT_BYTES {
        return Err(invalid("adduct text limit exceeded"));
    }
    Ok(())
}
fn finite(value: f64, what: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("nonfinite {what}")))
    }
}
fn consume(remaining: &mut usize, count: usize) -> Result<()> {
    *remaining = remaining
        .checked_sub(count)
        .ok_or_else(|| invalid("adduct cumulative work limit exceeded"))?;
    Ok(())
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
