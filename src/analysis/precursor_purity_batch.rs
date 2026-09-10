// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Acquisition-order precursor lookup and scalar purity for MS2 records.

use super::{
    MAX_PURITY_PRECURSORS, PrecursorPurity, PurityScores, WorkBudget, compute_with_budget, invalid,
    validate_tolerance,
};
use crate::comparison::Tolerance;
use crate::{MSExperiment, Result};
use std::collections::BTreeMap;

impl PrecursorPurity {
    /// Compute scalar purity for every MS2 native ID, using only its first
    /// precursor and the referenced or most recent preceding MS1 scan. Missing
    /// referenced IDs fall back to the most recent preceding parent. This does
    /// not average the preceding and following scans (the source header is stale).
    ///
    /// Missing parents yield zero scores only with ignore_missing=true. Invalid
    /// or duplicate MS2 IDs and otherwise missing required inputs return errors
    /// instead of the source warning plus empty map. Empty experiments are valid.
    /// At most 100,000 scans and a shared 50-million-unit work allowance apply.
    pub fn compute_all(
        experiment: &MSExperiment,
        tolerance: Tolerance,
        ignore_missing: bool,
    ) -> Result<BTreeMap<String, PurityScores>> {
        validate_tolerance(tolerance)?;
        if experiment.len() > MAX_PURITY_PRECURSORS {
            return Err(invalid("precursor purity experiment scan limit exceeded"));
        }
        if experiment.is_empty() {
            return Ok(BTreeMap::new());
        }
        if !ignore_missing && experiment.spectra[0].ms_level != 1 {
            return Err(invalid(
                "precursor purity experiment must start with an MS1 scan",
            ));
        }
        let mut budget = WorkBudget::new();
        let depth = (usize::BITS - experiment.len().leading_zeros()) as usize + 1;
        let mut by_id = BTreeMap::new();
        let mut latest = BTreeMap::new();
        let mut results = BTreeMap::new();
        for (index, spectrum) in experiment.spectra.iter().enumerate() {
            budget.consume(
                spectrum
                    .native_id
                    .len()
                    .saturating_add(1)
                    .saturating_mul(depth),
            )?;
            if spectrum.ms_level == 0 {
                return Err(invalid("MS level must be positive"));
            }
            if spectrum.ms_level == 2 {
                if spectrum.native_id.is_empty() || spectrum.native_id.chars().any(char::is_control)
                {
                    return Err(invalid(
                        "MS2 purity requires a nonempty control-free native ID",
                    ));
                }
                if results.contains_key(&spectrum.native_id) {
                    return Err(invalid(
                        "duplicate MS2 native ID in precursor purity experiment",
                    ));
                }
                let referenced = if let Some(reference) = spectrum
                    .precursors
                    .first()
                    .and_then(|p| p.spectrum_reference.as_deref())
                {
                    budget.consume(reference.len().saturating_add(1).saturating_mul(depth))?;
                    by_id.get(&(1, reference)).copied()
                } else {
                    None
                };
                let parent = referenced.or_else(|| latest.get(&1).copied());
                let score = if let Some(parent) = parent {
                    let precursor = spectrum
                        .precursors
                        .first()
                        .ok_or_else(|| invalid("MS2 spectrum has no precursor"))?;
                    compute_with_budget(
                        &experiment.spectra[parent],
                        precursor,
                        tolerance,
                        &mut budget,
                    )?
                } else if ignore_missing {
                    PurityScores::default()
                } else {
                    return Err(invalid("MS2 spectrum has no preceding parent spectrum"));
                };
                results.insert(spectrum.native_id.clone(), score);
            }
            latest.insert(spectrum.ms_level, index);
            by_id.insert((spectrum.ms_level, spectrum.native_id.as_str()), index);
        }
        Ok(results)
    }
}
