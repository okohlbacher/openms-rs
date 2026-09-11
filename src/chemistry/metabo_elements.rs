// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Narrow owned element-isotope projection for FeatureFindingMetabo.
use super::{EmpiricalFormula, element};
use crate::{Error, Result};
use std::mem::size_of;

pub(crate) fn isotope_masses(
    text: &str,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<Vec<Vec<f64>>> {
    let fail = || Error::InvalidValue("metabo element projection exceeds resource limits".into());
    // Formula input scan, ordered-map key comparisons and sparse nodes. The
    // alphabet is static, so even repeated source element names stay bounded.
    let count = text.len();
    *work = work
        .checked_sub(count.checked_mul(128).ok_or_else(fail)?)
        .ok_or_else(fail)?;
    *bytes = bytes
        .checked_sub(
            count
                .checked_mul(512)
                .and_then(|n| n.checked_add(512))
                .ok_or_else(fail)?,
        )
        .ok_or_else(fail)?;
    let formula = EmpiricalFormula::parse(text)?;
    let mut result = Vec::new();
    *bytes = bytes
        .checked_sub(
            formula
                .atoms
                .len()
                .checked_mul(size_of::<Vec<f64>>())
                .ok_or_else(fail)?,
        )
        .ok_or_else(fail)?;
    result
        .try_reserve_exact(formula.atoms.len())
        .map_err(|_| fail())?;
    for atom in formula.atoms.keys() {
        let isotopes = element(atom.symbol)
            .expect("parsed atom has a table")
            .isotopes();
        *work = work.checked_sub(isotopes.len()).ok_or_else(fail)?;
        let n = if atom.isotope.is_some() {
            1
        } else {
            isotopes.len()
        };
        *bytes = bytes
            .checked_sub(n.checked_mul(size_of::<f64>()).ok_or_else(fail)?)
            .ok_or_else(fail)?;
        let mut masses = Vec::new();
        masses.try_reserve_exact(n).map_err(|_| fail())?;
        for isotope in isotopes {
            if atom
                .isotope
                .is_none_or(|number| number == isotope.mass_number)
            {
                masses.push(isotope.mass);
            }
        }
        result.push(masses);
    }
    Ok(result)
}
